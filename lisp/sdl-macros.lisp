;;;; sdl-macros.lisp - SDL macros for simulation definition
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :sdl)

;;; ============================================================
;;; Global State
;;; ============================================================

(defvar *current-simulation* nil
  "The current simulation being defined")

(defun get-current-simulation ()
  "Get the current simulation object."
  *current-simulation*)

;;; ============================================================
;;; Simulation Structure
;;; ============================================================

(defstruct simulation-data
  "Internal structure holding simulation data"
  (name "unnamed" :type string)
  (version 1 :type integer)
  (dimension 2 :type integer)
  (geometry-shapes nil :type list)
  (mesh-spec nil)
  (mesh nil)
  (boundaries nil :type list)
  (physics nil)
  (discretization nil)
  (solver nil)
  (output nil))

;;; Accessors for C++ loader
(defun get-name (obj)
  "Get name from simulation-data or boundary-condition."
  (typecase obj
    (simulation-data (simulation-data-name obj))
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-name obj))
    (t nil)))

(defun get-version (sim)
  (simulation-data-version sim))

(defun get-dimension (sim)
  (simulation-data-dimension sim))

(defun get-mesh (sim)
  (simulation-data-mesh-spec sim))

(defun get-boundaries (sim)
  (simulation-data-boundaries sim))

(defun get-physics (sim)
  (simulation-data-physics sim))

(defun get-discretization (sim)
  (simulation-data-discretization sim))

(defun get-solver (sim)
  (simulation-data-solver sim))

(defun get-output (sim)
  (simulation-data-output sim))

;;; Generic property accessors for C++ bridge
;;; These handle different object types (mesh-spec, physics-data, etc.)

(defun get-type (obj)
  "Get type property from various objects."
  (typecase obj
    (physics-data (physics-data-type obj))
    (streamvorti.mesh:mesh-spec (streamvorti.mesh:mesh-spec-type obj))
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-type obj))
    (t nil)))

(defun get-path (obj)
  "Get path property from mesh-spec."
  (typecase obj
    (streamvorti.mesh:mesh-spec (streamvorti.mesh:mesh-spec-path obj))
    (t nil)))

(defun get-divisions (obj)
  "Get divisions property from mesh-spec."
  (typecase obj
    (streamvorti.mesh:mesh-spec (streamvorti.mesh:mesh-spec-divisions obj))
    (t nil)))

(defun get-element-type (obj)
  "Get element-type property from mesh-spec."
  (typecase obj
    (streamvorti.mesh:mesh-spec (streamvorti.mesh:mesh-spec-element-type obj))
    (t nil)))

(defun get-sizes (obj)
  "Get sizes property from mesh-spec."
  (typecase obj
    (streamvorti.mesh:mesh-spec (streamvorti.mesh:mesh-spec-sizes obj))
    (t nil)))

(defun get-size-x (obj)
  "Get size-x from mesh-spec (first element of sizes)."
  (typecase obj
    (streamvorti.mesh:mesh-spec
     (first (streamvorti.mesh:mesh-spec-sizes obj)))
    (t nil)))

(defun get-size-y (obj)
  "Get size-y from mesh-spec (second element of sizes)."
  (typecase obj
    (streamvorti.mesh:mesh-spec
     (second (streamvorti.mesh:mesh-spec-sizes obj)))
    (t nil)))

(defun get-size-z (obj)
  "Get size-z from mesh-spec (third element of sizes)."
  (typecase obj
    (streamvorti.mesh:mesh-spec
     (third (streamvorti.mesh:mesh-spec-sizes obj)))
    (t nil)))

(defun get-formulation (obj)
  "Get formulation from physics-data."
  (typecase obj
    (physics-data (physics-data-formulation obj))
    (t nil)))

(defun get-reynolds (obj)
  "Get Reynolds number from physics-data."
  (typecase obj
    (physics-data (physics-data-reynolds obj))
    (t nil)))

(defun get-density (obj)
  "Get density from physics-data."
  (typecase obj
    (physics-data (physics-data-density obj))
    (t nil)))

(defun get-viscosity (obj)
  "Get viscosity from physics-data."
  (typecase obj
    (physics-data (physics-data-viscosity obj))
    (t nil)))

(defun get-method (obj)
  "Get method from discretization-data."
  (typecase obj
    (discretization-data (discretization-data-method obj))
    (t nil)))

(defun get-num-neighbors (obj)
  "Get num-neighbors from discretization-data."
  (typecase obj
    (discretization-data (discretization-data-num-neighbors obj))
    (t nil)))

(defun get-cutoff-radius (obj)
  "Get cutoff-radius from discretization-data."
  (typecase obj
    (discretization-data (discretization-data-cutoff-radius obj))
    (t nil)))

(defun get-support-radius (obj)
  "Get support-radius from discretization-data."
  (typecase obj
    (discretization-data (discretization-data-support-radius obj))
    (t nil)))

(defun get-timestepping (obj)
  "Get timestepping from solver-data."
  (typecase obj
    (solver-data (solver-data-timestepping obj))
    (t nil)))

(defun get-dt (obj)
  "Get dt from solver-data."
  (typecase obj
    (solver-data (solver-data-dt obj))
    (t nil)))

(defun get-end-time (obj)
  "Get end-time from solver-data."
  (typecase obj
    (solver-data (solver-data-end-time obj))
    (t nil)))

(defun get-tolerance (obj)
  "Get tolerance from solver-data."
  (typecase obj
    (solver-data (solver-data-tolerance obj))
    (t nil)))

(defun get-max-iterations (obj)
  "Get max-iterations from solver-data."
  (typecase obj
    (solver-data (solver-data-max-iterations obj))
    (t nil)))

;;; Boundary condition accessors
(defun get-attribute (obj)
  "Get attribute from boundary-condition."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-attribute obj))
    (t nil)))

(defun get-function (obj)
  "Get function from boundary-condition or velocity plist."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-function obj))
    (cons (getf obj :function))  ; plist case
    (t nil)))

(defun get-u-function (obj)
  "Get u-function from boundary-condition or velocity plist."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-u-function obj))
    (cons (getf obj :u-function))
    (t nil)))

(defun get-v-function (obj)
  "Get v-function from boundary-condition or velocity plist."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-v-function obj))
    (cons (getf obj :v-function))
    (t nil)))

(defun get-w-function (obj)
  "Get w-function from boundary-condition or velocity plist."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-w-function obj))
    (cons (getf obj :w-function))
    (t nil)))

(defun get-predicate (obj)
  "Get predicate from boundary-condition."
  (typecase obj
    (streamvorti.boundaries:boundary-condition
     (streamvorti.boundaries:boundary-predicate obj))
    (t nil)))

(defun get-predicate-axis (obj)
  "Get the axis ('x', 'y', or 'z') from a simple coordinate predicate.
   Returns the axis string or empty string if not a simple predicate."
  (let ((pred (get-predicate obj)))
    (when pred
      (let ((axis (streamvorti.boundaries:predicate-axis pred)))
        (case axis
          (:x "x")
          (:y "y")
          (:z "z")
          (t ""))))))

(defun get-predicate-value (obj)
  "Get the coordinate value from a simple coordinate predicate.
   Returns the value or 0.0 if not determinable."
  (let ((pred (get-predicate obj)))
    (if pred
        (streamvorti.boundaries:predicate-value pred)
        0.0d0)))

(defun get-predicate-tolerance (obj)
  "Get the tolerance from a simple coordinate predicate."
  (let ((pred (get-predicate obj)))
    (if pred
        (streamvorti.boundaries:predicate-tolerance pred)
        1.0d-10)))

(defun evaluate-predicate (bc x y &optional (z 0.0d0))
  "Evaluate whether point (x,y,z) satisfies the boundary condition's predicate.
   Returns T if the predicate matches, NIL otherwise."
  (let ((pred (get-predicate bc)))
    (when pred
      (funcall (streamvorti.boundaries:predicate-test-function pred)
               x y z))))

;;; ============================================================
;;; Physics Structure
;;; ============================================================

(defstruct physics-data
  (type :incompressible-navier-stokes)
  (formulation :stream-vorticity)
  (reynolds 100.0d0 :type double-float)
  (density 1.0d0 :type double-float)
  (viscosity 0.01d0 :type double-float))

;;; ============================================================
;;; Discretization Structure
;;; ============================================================

(defstruct discretization-data
  (method :dcpse)
  (num-neighbors 25 :type integer)
  (cutoff-radius 30.0d0 :type double-float)
  (support-radius 5.0d0 :type double-float))

;;; ============================================================
;;; Solver Structure
;;; ============================================================

(defstruct solver-data
  (timestepping :explicit-euler)
  (dt 0.001d0 :type double-float)
  (end-time 1.0d0 :type double-float)
  (tolerance 1.0d-6 :type double-float)
  (max-iterations 10000 :type integer))

;;; ============================================================
;;; Output Structure
;;; ============================================================

(defstruct output-data
  (format :vtk)
  (interval 0.1d0 :type double-float)
  (directory "results/" :type string)
  (fields '(vorticity streamfunction velocity) :type list))

;;; Output data accessors (must be after defstruct)
(defun get-format (obj)
  "Get format from output-data."
  (typecase obj
    (output-data (output-data-format obj))
    (t nil)))

(defun get-interval (obj)
  "Get interval from output-data."
  (typecase obj
    (output-data (output-data-interval obj))
    (t nil)))

(defun get-directory (obj)
  "Get directory from output-data."
  (typecase obj
    (output-data (output-data-directory obj))
    (t nil)))

(defun get-fields (obj)
  "Get fields from output-data."
  (typecase obj
    (output-data (output-data-fields obj))
    (t nil)))

;;; ============================================================
;;; Main Simulation Macro
;;; ============================================================

(defmacro simulation (&rest args)
  "Define a simulation.

   Usage:
   (simulation
     :name \"cavity-flow\"
     :version 1
     :dimension 2

     (geometry ...)
     (mesh ...)
     (boundaries ...)
     (physics ...)
     (discretization ...)
     (solver ...)
     (output ...))"
  (let ((sim (gensym "SIM"))
        (name "unnamed")
        (version 1)
        (dimension 2)
        (body-clauses nil))
    ;; Parse args: extract keyword options and collect body clauses
    (loop with rest = args
          while rest
          for item = (car rest)
          do (cond
               ;; Keyword argument
               ((eq item :name)
                (setf name (cadr rest))
                (setf rest (cddr rest)))
               ((eq item :version)
                (setf version (cadr rest))
                (setf rest (cddr rest)))
               ((eq item :dimension)
                (setf dimension (cadr rest))
                (setf rest (cddr rest)))
               ;; Body clause (list starting with symbol)
               ((and (listp item) (symbolp (car item)))
                (push item body-clauses)
                (setf rest (cdr rest)))
               ;; Skip anything else
               (t
                (setf rest (cdr rest)))))
    (setf body-clauses (nreverse body-clauses))
    `(let ((,sim (make-simulation-data
                  :name ,name
                  :version ,version
                  :dimension ,dimension)))
       (setf *current-simulation* ,sim)
       ,@(mapcar (lambda (clause)
                   `(process-clause ,sim ',clause))
                 body-clauses)
       ,sim)))

(defgeneric process-clause (sim clause)
  (:documentation "Process a simulation clause"))

(defmethod process-clause (sim (clause list))
  "Process a clause based on its first element."
  (case (first clause)
    (geometry (process-geometry-clause sim (rest clause)))
    (mesh (process-mesh-clause sim (rest clause)))
    (boundaries (process-boundaries-clause sim (rest clause)))
    (physics (process-physics-clause sim (rest clause)))
    (discretization (process-discretization-clause sim (rest clause)))
    (solver (process-solver-clause sim (rest clause)))
    (output (process-output-clause sim (rest clause)))
    (otherwise
     (warn "Unknown clause type: ~A" (first clause)))))

;;; ============================================================
;;; Clause Processors
;;; ============================================================

(defun process-geometry-clause (sim body)
  "Process geometry definitions."
  ;; Execute body forms and collect defined shapes
  (dolist (form body)
    (when (and (listp form) (eq (first form) 'defparameter))
      (eval form)
      (let ((shape (symbol-value (second form))))
        (when (typep shape 'streamvorti.geometry:shape)
          (push shape (simulation-data-geometry-shapes sim)))))))

(defun process-mesh-clause (sim body)
  "Process mesh specification."
  (let ((mesh-spec nil))
    (dolist (form body)
      (cond
        ;; (generate shape ...)
        ((and (listp form) (eq (first form) 'generate))
         (let* ((geometry (eval (second form)))
                (options (cddr form))
                (type (or (getf options :type) :quad))
                ;; Evaluate divisions since it may be quoted '(40 40)
                (divisions-raw (or (getf options :divisions) '(10 10)))
                (divisions (if (and (listp divisions-raw)
                                    (eq (first divisions-raw) 'quote))
                               (second divisions-raw)
                               divisions-raw)))
           (setf mesh-spec
                 (streamvorti.mesh:make-mesh-spec-generate
                  geometry :type type :divisions divisions))))
        ;; (load path ...)
        ((and (listp form) (eq (first form) 'load-mesh))
         (let* ((path (second form))
                (format (or (getf (cddr form) :format) :auto)))
           (setf mesh-spec
                 (streamvorti.mesh:make-mesh-spec-load path :format format))))))
    (setf (simulation-data-mesh-spec sim) mesh-spec)))

(defun process-boundaries-clause (sim body)
  "Process boundary conditions."
  (let ((boundaries nil))
    (dolist (form body)
      (cond
        ;; (defun ...) - define BC function
        ((and (listp form) (eq (first form) 'defun))
         (eval form))
        ;; (region ...) - define BC region
        ((and (listp form) (eq (first form) 'region))
         (push (eval form) boundaries))))
    (setf (simulation-data-boundaries sim) (nreverse boundaries))))

(defun process-physics-clause (sim body)
  "Process physics parameters."
  (let ((physics (make-physics-data)))
    (loop for (key value) on body by #'cddr do
      (case key
        (:type (setf (physics-data-type physics) value))
        (:formulation (setf (physics-data-formulation physics) value))
        (:reynolds (setf (physics-data-reynolds physics) (coerce value 'double-float)))
        (:density (setf (physics-data-density physics) (coerce value 'double-float)))
        (:viscosity (setf (physics-data-viscosity physics) (coerce value 'double-float)))))
    (setf (simulation-data-physics sim) physics)))

(defun process-discretization-clause (sim body)
  "Process discretization parameters."
  (let ((disc (make-discretization-data)))
    (loop for (key value) on body by #'cddr do
      (case key
        (:method (setf (discretization-data-method disc) value))
        (:num-neighbors (setf (discretization-data-num-neighbors disc) value))
        (:cutoff-radius (setf (discretization-data-cutoff-radius disc)
                              (coerce value 'double-float)))
        (:support-radius (setf (discretization-data-support-radius disc)
                               (coerce value 'double-float)))))
    (setf (simulation-data-discretization sim) disc)))

(defun process-solver-clause (sim body)
  "Process solver parameters."
  (let ((solver (make-solver-data)))
    (loop for (key value) on body by #'cddr do
      (case key
        (:timestepping (setf (solver-data-timestepping solver) value))
        (:dt (setf (solver-data-dt solver) (coerce value 'double-float)))
        (:end-time (setf (solver-data-end-time solver) (coerce value 'double-float)))
        (:tolerance (setf (solver-data-tolerance solver) (coerce value 'double-float)))
        (:max-iterations (setf (solver-data-max-iterations solver) value))))
    (setf (simulation-data-solver sim) solver)))

(defun process-output-clause (sim body)
  "Process output configuration."
  (let ((output (make-output-data)))
    (loop for (key value) on body by #'cddr do
      (case key
        (:format (setf (output-data-format output) value))
        (:interval (setf (output-data-interval output) (coerce value 'double-float)))
        (:directory (setf (output-data-directory output) value))
        (:fields (setf (output-data-fields output) value))))
    (setf (simulation-data-output sim) output)))

;;; ============================================================
;;; Convenience Macros for SDL
;;; ============================================================

(defmacro geometry (&body body)
  "Define geometry shapes.
   This macro is used within simulation definition."
  `(progn ,@body))

(defmacro mesh (&body body)
  "Define mesh generation/loading.
   This macro is used within simulation definition."
  `(progn ,@body))

(defmacro boundaries (&body body)
  "Define boundary conditions.
   This macro is used within simulation definition."
  `(progn ,@body))

(defmacro physics (&body options)
  "Define physics parameters.
   This macro is used within simulation definition."
  (declare (ignore options))
  nil)

(defmacro discretization (&body options)
  "Define discretization parameters.
   This macro is used within simulation definition."
  (declare (ignore options))
  nil)

(defmacro solver (&body options)
  "Define solver parameters.
   This macro is used within simulation definition."
  (declare (ignore options))
  nil)

(defmacro output (&body options)
  "Define output configuration.
   This macro is used within simulation definition."
  (declare (ignore options))
  nil)

(defmacro generate (geometry &rest options)
  "Generate a mesh from geometry."
  `(streamvorti.mesh:generate-mesh ,geometry ,@options))

;;; ============================================================
;;; Simulation Loading
;;; ============================================================

(defun load-simulation (path)
  "Load a simulation definition from file."
  (setf *current-simulation* nil)
  (load path)
  (unless *current-simulation*
    (error "No simulation defined in ~A" path))
  *current-simulation*)
