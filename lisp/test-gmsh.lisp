;;;; test-gmsh.lisp - Test gmsh-cl integration
;;;;
;;;; Requires gmsh-cl to be loaded (STREAMVORTI_WITH_GMSH).
;;;; Run with:
;;;;   sbcl --load lisp/gmsh-init.lisp --load lisp/test-gmsh.lisp

(defpackage :streamvorti.gmsh.tests
  (:use :cl)
  (:export #:run-all-tests
           #:*pass-count*
           #:*fail-count*))

(in-package :streamvorti.gmsh.tests)

;;; ============================================================
;;; Test Framework (minimal, no dependencies)
;;; ============================================================

(defvar *test-count* 0)
(defvar *pass-count* 0)
(defvar *fail-count* 0)

(defmacro define-test (name &body body)
  "Define a test case."
  `(defun ,name ()
     (format t "~%Testing ~A...~%" ',name)
     (handler-case
         (progn ,@body t)
       (error (e)
         (incf *fail-count*)
         (format t "  ERROR: ~A~%" e)
         nil))))

(defun check (condition &optional (message ""))
  "Assert that CONDITION is true."
  (incf *test-count*)
  (if condition
      (incf *pass-count*)
      (progn
        (incf *fail-count*)
        (format t "  FAIL: ~A~%" message))))

(defun check-not-null (obj &optional (message ""))
  "Assert that OBJ is not NIL."
  (check (not (null obj)) message))

;;; ============================================================
;;; Tests
;;; ============================================================

(define-test test-gmsh-simple-rectangle
  "Generate a simple rectangle mesh with gmsh-cl."
  (let ((path (format nil "/tmp/streamvorti-test-~A.msh"
                      (get-universal-time))))
    (unwind-protect
        (progn
          (gmsh:with-gmsh ()
            (gmsh:with-model ("test-rect")
              (occ:rectangle 0 0 0 1 1)
              (occ:synchronize)
              (gmsh/mesh:generate :dim 2)
              (gmsh:write path)))
          (check (probe-file path) "Mesh file created")
          (check (> (with-open-file (s path) (file-length s)) 0)
                 "Mesh file is non-empty"))
      (when (probe-file path) (delete-file path)))))

(define-test test-gmsh-csg-difference
  "Generate a mesh with CSG boolean (rectangle minus circle)."
  (let ((path (format nil "/tmp/streamvorti-test-csg-~A.msh"
                      (get-universal-time))))
    (unwind-protect
        (progn
          (gmsh:with-gmsh ()
            (gmsh:with-model ("test-csg")
              (occ:rectangle 0 0 0 2 1 :tag 1)
              (occ:disk 0.5 0.5 0 0.1 0.1 :tag 2)
              (occ:cut (gmsh:surface-tags '(1))
                       (gmsh:surface-tags '(2)))
              (occ:synchronize)
              (gmsh/mesh:generate :dim 2)
              (gmsh:write path)))
          (check (probe-file path) "CSG mesh file created"))
      (when (probe-file path) (delete-file path)))))

(define-test test-gmsh-mesh-size-control
  "Verify mesh size option affects element count."
  (let ((path-coarse (format nil "/tmp/sv-test-coarse-~A.msh"
                             (get-universal-time)))
        (path-fine (format nil "/tmp/sv-test-fine-~A.msh"
                           (get-universal-time))))
    (unwind-protect
        (progn
          ;; Coarse mesh
          (gmsh:with-gmsh ()
            (gmsh:with-model ("coarse")
              (occ:rectangle 0 0 0 1 1)
              (occ:synchronize)
              (gmsh/option:set-number "Mesh.CharacteristicLengthMax" 0.5d0)
              (gmsh/mesh:generate :dim 2)
              (gmsh:write path-coarse)))
          ;; Fine mesh
          (gmsh:with-gmsh ()
            (gmsh:with-model ("fine")
              (occ:rectangle 0 0 0 1 1)
              (occ:synchronize)
              (gmsh/option:set-number "Mesh.CharacteristicLengthMax" 0.1d0)
              (gmsh/mesh:generate :dim 2)
              (gmsh:write path-fine)))
          ;; Fine should produce a larger file
          (let ((coarse-size (with-open-file (s path-coarse) (file-length s)))
                (fine-size (with-open-file (s path-fine) (file-length s))))
            (check (> fine-size coarse-size)
                   "Fine mesh file larger than coarse")))
      (when (probe-file path-coarse) (delete-file path-coarse))
      (when (probe-file path-fine) (delete-file path-fine)))))

(define-test test-gmsh-domain-in-simulation
  "Test (domain :gmsh ...) with physical groups and attribute predicates."
  (let ((sim (sdl:simulation "gmsh-sim-test" :dim 2
               (domain :gmsh
                 (occ:rectangle 0 0 0 1 1)
                 (occ:synchronize)
                 (let ((bnd (gmsh:get-boundary (gmsh:get-entities :dim 2)
                              :oriented nil)))
                   (gmsh:add-physical-group 1
                     (gmsh:tags-of (occ:get-closest-entities 0 0.5 0 bnd :n 1))
                     :tag 1 :name "left")
                   (gmsh:add-physical-group 1
                     (gmsh:tags-of (occ:get-closest-entities 1 0.5 0 bnd :n 1))
                     :tag 2 :name "right")
                   (gmsh:add-physical-group 1
                     (gmsh:tags-of (occ:get-closest-entities 0.5 0 0 bnd :n 1))
                     :tag 3 :name "bottom")
                   (gmsh:add-physical-group 1
                     (gmsh:tags-of (occ:get-closest-entities 0.5 1 0 bnd :n 1))
                     :tag 4 :name "top"))
                 (gmsh:add-physical-group 2
                   (gmsh:tags-of (gmsh:get-entities :dim 2)) :tag 1)
                 (gmsh/mesh:generate :dim 2))
               (boundaries
                 (left   (attribute 1))
                 (right  (attribute 2))
                 (bottom (attribute 3))
                 (top    (attribute 4)))
               (physics :navier-stokes :Re 100
                 (bc left   :no-slip)
                 (bc right  :no-slip)
                 (bc bottom :no-slip)
                 (bc top    :velocity (1 0)))
               (temporal :explicit-euler :dt 0.001 :end 0.01))))
    (check-not-null sim "Gmsh simulation created")
    (check-not-null (sdl:simulation-domain sim) "Domain defined")
    (let ((domain (sdl:simulation-domain sim)))
      (check (eq (sdl:domain-type domain) :file)
             "Domain type is :file")
      (check-not-null (sdl:domain-file domain)
                      "Domain has file path"))))

;;; ============================================================
;;; Run All Tests
;;; ============================================================

(defun run-all-tests ()
  "Run all gmsh-cl integration tests."
  (setf *test-count* 0
        *pass-count* 0
        *fail-count* 0)

  (format t "~%========================================~%")
  (format t "StreamVorti Gmsh Integration Tests~%")
  (format t "========================================~%")

  (test-gmsh-simple-rectangle)
  (test-gmsh-csg-difference)
  (test-gmsh-mesh-size-control)
  (test-gmsh-domain-in-simulation)

  (format t "~%========================================~%")
  (format t "Results: ~D passed, ~D failed, ~D total~%"
          *pass-count* *fail-count* *test-count*)
  (format t "========================================~%")

  (zerop *fail-count*))

;;; Run tests when loaded
(unless (run-all-tests)
  (format t "~%GMSH TESTS FAILED~%"))
