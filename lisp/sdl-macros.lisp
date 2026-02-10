;;;; sdl-macros.lisp - SDL v2 API implementation
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :sdl)

;; SBCL considers make-method a reserved CL name (method is a standard class).
;; Unlock the CL package to allow defining it.
#+sbcl (sb-ext:unlock-package :common-lisp)

;;; ============================================================
;;; Geometry Wrappers
;;; ============================================================

(defun box (min-corner max-corner)
  "Create a box geometry, dispatching by dimension.
   1D (single-element lists) -> interval
   2D (two-element lists)    -> rectangle
   3D (three-element lists)  -> box"
  (ecase (length min-corner)
    (1 (make-instance 'interval
                      :min-corner (coerce min-corner 'list)
                      :max-corner (coerce max-corner 'list)))
    (2 (make-rectangle min-corner max-corner))
    (3 (make-box min-corner max-corner))))

(defun ball (center radius)
  "Create a ball geometry, dispatching by dimension.
   2D (two-element center) -> circle
   3D (three-element center) -> sphere"
  (ecase (length center)
    (2 (make-circle center radius))
    (3 (make-sphere center radius))))

(defun geometry-dimension (geom)
  "Return the spatial dimension of a geometry."
  (shape-dimension geom))

(defun geometry-min (geom)
  "Return the minimum corner of a geometry as a list."
  (typecase geom
    (interval (interval-min-corner geom))
    (rectangle (rectangle-min-corner geom))
    (streamvorti.geometry:box (streamvorti.geometry:box-min-corner geom))))

(defun geometry-max (geom)
  "Return the maximum corner of a geometry as a list."
  (typecase geom
    (interval (interval-max-corner geom))
    (rectangle (rectangle-max-corner geom))
    (streamvorti.geometry:box (streamvorti.geometry:box-max-corner geom))))

(defun geometry-contains-p (geom point)
  "Test if a geometry contains a point."
  (shape-contains-p geom point))

(defun ball-center (ball)
  "Return the center of a ball geometry."
  (typecase ball
    (circle (circle-center ball))
    (sphere (sphere-center ball))))

(defun ball-radius (ball)
  "Return the radius of a ball geometry."
  (typecase ball
    (circle (circle-radius ball))
    (sphere (sphere-radius ball))))

;;; ============================================================
;;; Domain
;;; ============================================================

(defstruct (domain-data (:constructor %make-domain)
                        (:conc-name domain-))
  type dimension geometry n h file particle-count)

(defun domain (&rest args)
  "Create a domain.

   Usage:
   (domain geometry :mesh :n '(10 10))
   (domain geometry :points :n '(20 20))
   (domain geometry :points :h 0.05)
   (domain geometry :particles :n 1000)
   (domain :file \"path.mesh\")"
  (cond
    ;; File domain: (domain :file "path")
    ((eq (first args) :file)
     (%make-domain :type :file :file (second args)))
    ;; Geometry-based domain: (domain geom type &key ...)
    (t
     (let* ((geom (first args))
            (type (second args))
            (rest (cddr args))
            (n nil) (h nil) (particle-count nil))
       ;; Parse keyword args
       (loop with r = rest
             while r
             do (cond
                  ((eq (car r) :n)
                   (if (eq type :particles)
                       (setf particle-count (cadr r))
                       (setf n (cadr r)))
                   (setf r (cddr r)))
                  ((eq (car r) :h)
                   (setf h (cadr r))
                   (setf r (cddr r)))
                  (t (setf r (cdr r)))))
       (%make-domain
        :type type
        :dimension (shape-dimension geom)
        :geometry geom
        :n n
        :h (when h (coerce h 'double-float))
        :particle-count particle-count)))))

;;; ============================================================
;;; Predicates
;;; ============================================================

(defstruct (predicate-data (:constructor %make-predicate)
                           (:conc-name predicate-))
  test-fn)

(defun coord-index (sym)
  "Map a coordinate symbol to its index in a point list."
  (let ((name (string-upcase (symbol-name sym))))
    (cond
      ((string= name "X") 0)
      ((string= name "Y") 1)
      ((string= name "Z") 2)
      (t (error "Unknown coordinate: ~A" sym)))))

(defun compile-predicate-expr (expr)
  "Compile an S-expression predicate to a test function (point) -> boolean."
  (let ((op (first expr)))
    (cond
      ;; (= coord value)
      ((string= (symbol-name op) "=")
       (let ((idx (coord-index (second expr)))
             (val (coerce (third expr) 'double-float)))
         (lambda (point)
           (< (abs (- (nth idx point) val)) 1.0d-10))))
      ;; (> coord value)
      ((string= (symbol-name op) ">")
       (let ((idx (coord-index (second expr)))
             (val (coerce (third expr) 'double-float)))
         (lambda (point)
           (> (nth idx point) val))))
      ;; (< coord value)
      ((string= (symbol-name op) "<")
       (let ((idx (coord-index (second expr)))
             (val (coerce (third expr) 'double-float)))
         (lambda (point)
           (< (nth idx point) val))))
      ;; (>= coord value)
      ((string= (symbol-name op) ">=")
       (let ((idx (coord-index (second expr)))
             (val (coerce (third expr) 'double-float)))
         (lambda (point)
           (>= (nth idx point) val))))
      ;; (<= coord value)
      ((string= (symbol-name op) "<=")
       (let ((idx (coord-index (second expr)))
             (val (coerce (third expr) 'double-float)))
         (lambda (point)
           (<= (nth idx point) val))))
      ;; (and ...)
      ((string= (symbol-name op) "AND")
       (let ((sub-fns (mapcar #'compile-predicate-expr (rest expr))))
         (lambda (point)
           (every (lambda (fn) (funcall fn point)) sub-fns))))
      ;; (or ...)
      ((string= (symbol-name op) "OR")
       (let ((sub-fns (mapcar #'compile-predicate-expr (rest expr))))
         (lambda (point)
           (some (lambda (fn) (funcall fn point)) sub-fns))))
      ;; (not ...)
      ((string= (symbol-name op) "NOT")
       (let ((sub-fn (compile-predicate-expr (second expr))))
         (lambda (point)
           (not (funcall sub-fn point)))))
      (t (error "Unknown predicate operator: ~A" op)))))

(defun make-predicate (expr)
  "Create a predicate from an S-expression.

   Supported forms:
   (= x 0)    (= y 1)    (= z 0.5)
   (> x 0)    (< x 1)    (>= x 0)    (<= x 1)
   (and ...)  (or ...)   (not ...)"
  (%make-predicate :test-fn (compile-predicate-expr expr)))

(defun predicate-matches-p (pred point)
  "Test if a point matches a predicate."
  (funcall (predicate-test-fn pred) point))

;;; ============================================================
;;; Boundary Selectors
;;; ============================================================

(defstruct (boundary-selector-data (:constructor %make-boundary-selector)
                                   (:conc-name selector-))
  type attribute)

(defun make-boundary-selector (spec)
  "Create a boundary selector from a specification.
   Example: (make-boundary-selector '(attribute 3))"
  (let ((type-name (string-upcase (symbol-name (first spec)))))
    (cond
      ((string= type-name "ATTRIBUTE")
       (%make-boundary-selector :type :attribute :attribute (second spec)))
      (t (error "Unknown selector type: ~A" (first spec))))))

;;; ============================================================
;;; Boundaries
;;; ============================================================

(defstruct (boundary-def (:constructor %make-boundary-def)
                         (:conc-name boundary-))
  name predicate)

(defstruct (boundary-set (:constructor %make-boundary-set))
  defs)

(defun make-boundaries (specs)
  "Create boundary definitions from a specification list.

   Example:
   (make-boundaries '((lid (= y 1))
                      (walls (or (= x 0) (= x 1)))))"
  (mapcar (lambda (spec)
            (%make-boundary-def
             :name (first spec)
             :predicate (make-predicate (second spec))))
          specs))

;;; ============================================================
;;; Subdomains
;;; ============================================================

(defstruct (subdomain-def (:constructor %make-subdomain-def)
                          (:conc-name subdomain-))
  name predicate)

(defstruct (subdomain-set (:constructor %make-subdomain-set))
  defs)

(defun make-subdomains (specs)
  "Create subdomain definitions from a specification list.

   Example:
   (make-subdomains '((fluid (< x 0.5))
                      (solid (>= x 0.5))))"
  (mapcar (lambda (spec)
            (%make-subdomain-def
             :name (first spec)
             :predicate (make-predicate (second spec))))
          specs))

;;; ============================================================
;;; Boundary Conditions
;;; ============================================================

(defstruct (bc-data (:constructor %make-bc)
                    (:conc-name bc-))
  type boundary value function h T-inf)

(defun make-bc (boundary type &rest args)
  "Create a boundary condition.

   Examples:
   (make-bc 'inlet :velocity '(1 0 0))
   (make-bc 'wall :no-slip)
   (make-bc 'outlet :pressure 0)
   (make-bc 'surface :convection :h 10 :T-inf 20)
   (make-bc 'inlet :velocity :fn (lambda (x y z) ...))"
  (let ((value nil) (fn nil) (h nil) (T-inf nil))
    ;; Parse args
    (loop with rest = args
          while rest
          do (cond
               ((eq (car rest) :fn)
                (setf fn (cadr rest))
                (setf rest (cddr rest)))
               ((eq (car rest) :h)
                (setf h (cadr rest))
                (setf rest (cddr rest)))
               ((eq (car rest) :T-inf)
                (setf T-inf (cadr rest))
                (setf rest (cddr rest)))
               ;; Positional value
               (t
                (setf value (car rest))
                (setf rest (cdr rest)))))
    (%make-bc
     :boundary boundary
     :type type
     :value value
     :function fn
     :h (when h (coerce h 'double-float))
     :T-inf (when T-inf (coerce T-inf 'double-float)))))

;;; ============================================================
;;; Physics
;;; ============================================================

(defstruct (physics-data (:constructor %make-physics)
                         (:conc-name physics-))
  name type (Re nil) (conductivity nil) (formulation nil) subdomain bcs)

(defun parse-bc-spec (spec)
  "Parse a boundary condition specification like (lid :velocity (1 0))."
  (let ((boundary (first spec))
        (type (second spec))
        (value (third spec)))
    (%make-bc :boundary boundary :type type :value value)))

(defun make-physics (&key type Re conductivity bcs)
  "Create a physics definition.

   Example:
   (make-physics :type :navier-stokes :Re 100
                 :bcs '((lid :velocity (1 0))
                        (walls :no-slip)))"
  (%make-physics
   :type type
   :Re (when Re (coerce Re 'double-float))
   :conductivity (when conductivity (coerce conductivity 'double-float))
   :bcs (when bcs (mapcar #'parse-bc-spec bcs))))

;;; ============================================================
;;; Methods (Discretization)
;;; ============================================================

(defstruct (method-data (:constructor %make-method)
                        (:conc-name method-))
  type neighbors support-radius order kernel h)

(defun make-method (type &key neighbors support-radius order kernel h)
  "Create a discretization method.

   Examples:
   (make-method :dcpse :neighbors 25 :support-radius 3.5)
   (make-method :fem :order 2)
   (make-method :sph :kernel :wendland :h 0.02)"
  (%make-method
   :type type
   :neighbors neighbors
   :support-radius (when support-radius (coerce support-radius 'double-float))
   :order order
   :kernel kernel
   :h (when h (coerce h 'double-float))))

;;; ============================================================
;;; Solvers
;;; ============================================================

(defstruct (solver-data (:constructor %make-solver)
                        (:conc-name solver-))
  method dt end tolerance)

(defun make-time-solver (&key method dt end)
  "Create a time integration solver.

   Example:
   (make-time-solver :method :bdf2 :dt 0.001 :end 10.0)"
  (%make-solver
   :method method
   :dt (when dt (coerce dt 'double-float))
   :end (when end (coerce end 'double-float))))

(defun make-steady-solver (&key method tolerance)
  "Create a steady-state solver.

   Example:
   (make-steady-solver :method :newton :tolerance 1e-8)"
  (%make-solver
   :method method
   :tolerance (when tolerance (coerce tolerance 'double-float))))

(defstruct (linear-solver-data (:constructor %make-linear-solver)
                               (:conc-name linear-solver-))
  method preconditioner)

(defun make-linear-solver (&key method preconditioner)
  "Create a linear solver configuration.

   Example:
   (make-linear-solver :method :gmres :preconditioner :ilu)"
  (%make-linear-solver
   :method method
   :preconditioner preconditioner))

;;; ============================================================
;;; Materials
;;; ============================================================

(defstruct (material-data (:constructor %make-material)
                          (:conc-name material-))
  name density youngs-modulus poissons-ratio viscosity conductivity)

(defun make-material (name &key density youngs-modulus poissons-ratio
                               viscosity conductivity)
  "Create a material definition.

   Examples:
   (make-material 'steel :density 7850 :youngs-modulus 200e9 :poissons-ratio 0.3)
   (make-material 'water :density 1000 :viscosity 0.001 :conductivity 0.6)"
  (%make-material
   :name name
   :density (when density (coerce density 'double-float))
   :youngs-modulus (when youngs-modulus (coerce youngs-modulus 'double-float))
   :poissons-ratio (when poissons-ratio (coerce poissons-ratio 'double-float))
   :viscosity (when viscosity (coerce viscosity 'double-float))
   :conductivity (when conductivity (coerce conductivity 'double-float))))

;;; ============================================================
;;; Output
;;; ============================================================

(defstruct (output-data (:constructor %make-output)
                        (:conc-name output-))
  format directory interval fields probes)

(defun make-output (&key format directory every fields probes)
  "Create output configuration.

   Example:
   (make-output :format :vtk :directory \"results/\" :every 0.1
                :fields '(velocity pressure))"
  (%make-output
   :format format
   :directory directory
   :interval (when every (coerce every 'double-float))
   :fields fields
   :probes probes))

;;; ============================================================
;;; Coupling (Multiphysics)
;;; ============================================================

(defstruct (coupling-data (:constructor %make-coupling)
                          (:conc-name coupling-))
  name physics type iterations interface)

(defun make-coupling (&key name physics type iterations interface)
  "Create a multiphysics coupling definition.

   Examples:
   (make-coupling :name :thermal-flow :physics '(thermal flow)
                  :type :sequential :iterations 10)
   (make-coupling :name :fsi :physics '(fluid structure)
                  :type :monolithic :interface 'moving-wall)"
  (%make-coupling
   :name name
   :physics physics
   :type type
   :iterations iterations
   :interface interface))

#+sbcl (sb-ext:lock-package :common-lisp)
