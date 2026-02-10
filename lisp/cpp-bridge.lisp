;;;; cpp-bridge.lisp - Bridge between C++ Loader and SDL v2 data structures
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.
;;;;
;;;; The C++ Loader::getProperty() method constructs Lisp function names
;;;; as "sdl:get-<key>" and calls them via Bridge::funcall().  This file
;;;; defines those accessor functions so the C++ code can extract
;;;; configuration from the SDL v2 data structures.

(in-package :sdl)

;;; ============================================================
;;; Merged Boundary Condition (combines boundary-def + bc-data)
;;; ============================================================

(defstruct (merged-bc (:conc-name merged-bc-))
  name attribute type-str
  predicate-axis predicate-value (predicate-tolerance 1e-10)
  function u-function v-function w-function)

;;; ============================================================
;;; Predicate analysis helpers
;;; ============================================================

(defun simple-equality-predicate-p (expr)
  "Return T if EXPR is a simple (= coord value) predicate."
  (and (listp expr)
       (= (length expr) 3)
       (symbolp (first expr))
       (string= (symbol-name (first expr)) "=")
       (symbolp (second expr))))

(defun extract-predicate-axis (expr)
  "Extract axis string from simple (= coord val) predicate, or nil."
  (when (simple-equality-predicate-p expr)
    (string-downcase (symbol-name (second expr)))))

(defun extract-predicate-value (expr)
  "Extract value from simple (= coord val) predicate, or nil."
  (when (simple-equality-predicate-p expr)
    (coerce (third expr) 'double-float)))

;;; ============================================================
;;; Boundary merger
;;; ============================================================

(defun merge-sim-boundaries (sim)
  "Create merged BC list combining boundary-defs with physics BCs.
   Returns a proper Lisp list of merged-bc structs."
  (let ((boundaries (sim-boundaries sim))
        (physics-list (sim-physics-list sim))
        (result nil)
        (attr 1))
    (when boundaries
      ;; Collect all BCs from all physics
      (let ((all-bcs nil))
        (dolist (phys physics-list)
          (dolist (bc (physics-bcs phys))
            (push bc all-bcs)))
        (setf all-bcs (nreverse all-bcs))
        ;; Match boundaries with BCs
        (dolist (bdef boundaries)
          (let* ((bname (boundary-name bdef))
                 (bname-str (string-upcase (symbol-name bname)))
                 (bc (find-if (lambda (bc)
                                (string= bname-str
                                         (string-upcase
                                          (symbol-name (bc-boundary bc)))))
                              all-bcs))
                 (pexpr (boundary-predicate-expr bdef))
                 (merged (make-merged-bc
                          :name (string-downcase (symbol-name bname))
                          :attribute attr
                          :type-str (if bc
                                        (string-downcase
                                         (symbol-name (bc-type bc)))
                                        "no-slip")
                          :predicate-axis (extract-predicate-axis pexpr)
                          :predicate-value (extract-predicate-value pexpr))))
            ;; Set up BC functions
            (when bc
              (cond
                ;; BC has a function (combined lambda returning list)
                ((bc-function bc)
                 (let ((f (bc-function bc)))
                   (setf (merged-bc-function merged) f)
                   ;; Create component wrappers
                   (setf (merged-bc-u-function merged)
                         (lambda (x y &optional z)
                           (first (funcall f x y z))))
                   (setf (merged-bc-v-function merged)
                         (lambda (x y &optional z)
                           (second (funcall f x y z))))))
                ;; BC has a constant value list like (1 0)
                ((and (bc-value bc) (listp (bc-value bc)))
                 (let ((u-val (coerce (first (bc-value bc)) 'double-float))
                       (v-val (coerce (second (bc-value bc)) 'double-float)))
                   (setf (merged-bc-u-function merged)
                         (lambda (x y &optional z)
                           (declare (ignore x y z))
                           u-val))
                   (setf (merged-bc-v-function merged)
                         (lambda (x y &optional z)
                           (declare (ignore x y z))
                           v-val))))))
            (push merged result)
            (incf attr)))))
    (nreverse result)))

;;; ============================================================
;;; Simulation accessors (sdl:get-* for sim-data)
;;; ============================================================

(defun get-name (obj)
  "Get name from simulation or merged-BC."
  (typecase obj
    (sim-data (sim-name obj))
    (merged-bc (merged-bc-name obj))
    (t nil)))

(defun get-version (obj)
  "Get version number (always 2 for SDL v2)."
  (declare (ignore obj))
  2)

(defun get-dimension (obj)
  "Get dimension from simulation or domain."
  (typecase obj
    (sim-data (sim-dim obj))
    (domain-data (domain-dimension obj))
    (t nil)))

(defun get-mesh (obj)
  "Get mesh/domain specification from simulation."
  (typecase obj
    (sim-data (sim-domain obj))
    (t nil)))

(defun get-boundaries (obj)
  "Get merged boundary conditions from simulation."
  (typecase obj
    (sim-data (merge-sim-boundaries obj))
    (t nil)))

(defun get-physics (obj)
  "Get primary physics from simulation."
  (typecase obj
    (sim-data (first (sim-physics-list obj)))
    (t nil)))

(defun get-discretization (obj)
  "Get spatial discretization from simulation."
  (typecase obj
    (sim-data (sim-spatial obj))
    (t nil)))

(defun get-solver (obj)
  "Get temporal/solver config from simulation."
  (typecase obj
    (sim-data (sim-temporal obj))
    (t nil)))

(defun get-output (obj)
  "Get output config from simulation."
  (typecase obj
    (sim-data (sim-output-config obj))
    (t nil)))

;;; ============================================================
;;; Domain/mesh accessors
;;; ============================================================

(defun get-type (obj)
  "Get type from domain, physics, or boundary condition."
  (typecase obj
    (domain-data
     (case (domain-type obj)
       ((:mesh :points) "generated")
       (:file "loaded")
       (t (string-downcase (symbol-name (domain-type obj))))))
    (physics-data
     (string-downcase (symbol-name (physics-type obj))))
    (merged-bc
     (merged-bc-type-str obj))
    (t nil)))

(defun get-divisions (obj)
  "Get mesh divisions from domain."
  (typecase obj
    (domain-data (domain-n obj))
    (t nil)))

(defun get-size-x (obj)
  "Get x-size from domain geometry."
  (typecase obj
    (domain-data
     (let ((geom (domain-geometry obj)))
       (when geom
         (let ((max-pt (geometry-max geom))
               (min-pt (geometry-min geom)))
           (coerce (- (first max-pt) (first min-pt)) 'double-float)))))
    (t nil)))

(defun get-size-y (obj)
  "Get y-size from domain geometry."
  (typecase obj
    (domain-data
     (let ((geom (domain-geometry obj)))
       (when geom
         (let ((max-pt (geometry-max geom))
               (min-pt (geometry-min geom)))
           (when (> (length max-pt) 1)
             (coerce (- (second max-pt) (second min-pt))
                     'double-float))))))
    (t nil)))

(defun get-size-z (obj)
  "Get z-size from domain geometry."
  (typecase obj
    (domain-data
     (let ((geom (domain-geometry obj)))
       (when geom
         (let ((max-pt (geometry-max geom))
               (min-pt (geometry-min geom)))
           (when (> (length max-pt) 2)
             (coerce (- (third max-pt) (third min-pt))
                     'double-float))))))
    (t nil)))

(defun get-element-type (obj)
  "Get element type string from domain."
  (typecase obj
    (domain-data
     (if (and (domain-dimension obj)
              (= (domain-dimension obj) 3))
         "hex"
         "quad"))
    (t nil)))

(defun get-path (obj)
  "Get file path from domain."
  (typecase obj
    (domain-data (domain-file obj))
    (t nil)))

;;; ============================================================
;;; Merged-BC accessors
;;; ============================================================

(defun get-attribute (obj)
  "Get boundary attribute number."
  (typecase obj
    (merged-bc (merged-bc-attribute obj))
    (t nil)))

(defun get-predicate-axis (obj)
  "Get predicate axis string."
  (typecase obj
    (merged-bc (merged-bc-predicate-axis obj))
    (t nil)))

(defun get-predicate-value (obj)
  "Get predicate coordinate value."
  (typecase obj
    (merged-bc (merged-bc-predicate-value obj))
    (t nil)))

(defun get-predicate-tolerance (obj)
  "Get predicate tolerance."
  (typecase obj
    (merged-bc (merged-bc-predicate-tolerance obj))
    (t nil)))

(defun get-function (obj)
  "Get combined BC function."
  (typecase obj
    (merged-bc (merged-bc-function obj))
    (t nil)))

(defun get-u-function (obj)
  "Get u-component BC function."
  (typecase obj
    (merged-bc (merged-bc-u-function obj))
    (t nil)))

(defun get-v-function (obj)
  "Get v-component BC function."
  (typecase obj
    (merged-bc (merged-bc-v-function obj))
    (t nil)))

(defun get-w-function (obj)
  "Get w-component BC function."
  (typecase obj
    (merged-bc (merged-bc-w-function obj))
    (t nil)))

;;; ============================================================
;;; Physics accessors
;;; ============================================================

(defun get-formulation (obj)
  "Get physics formulation string."
  (typecase obj
    (physics-data
     (if (physics-formulation obj)
         (string-downcase (symbol-name (physics-formulation obj)))
         "stream-vorticity"))
    (t nil)))

(defun get-reynolds (obj)
  "Get Reynolds number."
  (typecase obj
    (physics-data (physics-Re obj))
    (t nil)))

(defun get-density (obj)
  "Get density (not stored in current structs)."
  (declare (ignore obj))
  nil)

(defun get-viscosity (obj)
  "Get viscosity (not stored in current structs)."
  (declare (ignore obj))
  nil)

;;; ============================================================
;;; Spatial/discretization accessors
;;; ============================================================

(defun get-num-neighbors (obj)
  "Get number of neighbors from spatial config."
  (typecase obj
    (spatial-data (spatial-neighbors obj))
    (t nil)))

(defun get-cutoff-radius (obj)
  "Get cutoff radius (not stored in current structs)."
  (declare (ignore obj))
  nil)

(defun get-support-radius (obj)
  "Get support radius from spatial config."
  (typecase obj
    (spatial-data (spatial-support-radius obj))
    (t nil)))

;;; ============================================================
;;; Temporal/solver accessors
;;; ============================================================

(defun get-timestepping (obj)
  "Get time-stepping method string."
  (typecase obj
    (temporal-data
     (if (temporal-method obj)
         (string-downcase (symbol-name (temporal-method obj)))
         "explicit-euler"))
    (t nil)))

(defun get-dt (obj)
  "Get time step size."
  (typecase obj
    (temporal-data (temporal-dt obj))
    (t nil)))

(defun get-end-time (obj)
  "Get simulation end time."
  (typecase obj
    (temporal-data (temporal-end obj))
    (t nil)))

(defun get-tolerance (obj)
  "Get solver tolerance."
  (typecase obj
    (temporal-data (temporal-tolerance obj))
    (t nil)))

(defun get-max-iterations (obj)
  "Get maximum iterations (not stored in current structs)."
  (declare (ignore obj))
  nil)

;;; ============================================================
;;; Output accessors
;;; ============================================================

(defun get-format (obj)
  "Get output format string."
  (typecase obj
    (output-data
     (if (output-format obj)
         (string-downcase (symbol-name (output-format obj)))
         nil))
    (t nil)))

(defun get-interval (obj)
  "Get output interval."
  (typecase obj
    (output-data (output-interval obj))
    (t nil)))

(defun get-directory (obj)
  "Get output directory."
  (typecase obj
    (output-data (output-directory obj))
    (t nil)))

(defun get-fields (obj)
  "Get output fields list."
  (typecase obj
    (output-data (output-fields obj))
    (t nil)))
