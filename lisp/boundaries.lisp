;;;; boundaries.lisp - Boundary condition handling
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :streamvorti.boundaries)

;;; ============================================================
;;; Boundary Predicates
;;; ============================================================

(defclass boundary-predicate ()
  ((test-function :initarg :test-function
                  :accessor predicate-test-function
                  :documentation "Function (x y z) -> boolean"))
  (:documentation "Predicate for selecting boundary elements"))

(defun x-equals (value &optional (tolerance 1e-10))
  "Create predicate for x = value."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (declare (ignore y z))
                                  (< (abs (- x value)) tolerance))))

(defun y-equals (value &optional (tolerance 1e-10))
  "Create predicate for y = value."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (declare (ignore x z))
                                  (< (abs (- y value)) tolerance))))

(defun z-equals (value &optional (tolerance 1e-10))
  "Create predicate for z = value."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (declare (ignore x y))
                                  (< (abs (- z value)) tolerance))))

(defun and-predicate (&rest predicates)
  "Create predicate that is true when all sub-predicates are true."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (every (lambda (p)
                                           (funcall (predicate-test-function p) x y z))
                                         predicates))))

(defun or-predicate (&rest predicates)
  "Create predicate that is true when any sub-predicate is true."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (some (lambda (p)
                                          (funcall (predicate-test-function p) x y z))
                                        predicates))))

(defun not-predicate (predicate)
  "Create predicate that negates another predicate."
  (make-instance 'boundary-predicate
                 :test-function (lambda (x y z)
                                  (not (funcall (predicate-test-function predicate)
                                                x y z)))))

;;; Helper macro for where clause
(defmacro where (expr)
  "Create a boundary predicate from an expression.

   Supported forms:
   (where (= x 0.0))   -> x-equals
   (where (= y 1.0))   -> y-equals
   (where (= z 0.5))   -> z-equals
   (where (and ...))   -> and-predicate
   (where (or ...))    -> or-predicate"
  (labels ((transform (e)
             (cond
               ;; (= x value) or (= y value) or (= z value)
               ((and (listp e) (eq (first e) '=))
                (let ((var (second e))
                      (val (third e)))
                  ;; Compare symbol names to handle cross-package usage
                  (let ((var-name (symbol-name var)))
                    (cond
                      ((string= var-name "X") `(x-equals ,val))
                      ((string= var-name "Y") `(y-equals ,val))
                      ((string= var-name "Z") `(z-equals ,val))
                      (t (error "Unknown variable in where: ~A" var))))))
               ;; (< x value) etc - not yet implemented
               ((and (listp e) (member (first e) '(< > <= >=)))
                (error "Comparison ~A not yet implemented in where" (first e)))
               ;; (and ...)
               ((and (listp e) (eq (first e) 'and))
                `(and-predicate ,@(mapcar #'transform (rest e))))
               ;; (or ...)
               ((and (listp e) (eq (first e) 'or))
                `(or-predicate ,@(mapcar #'transform (rest e))))
               ;; (not ...)
               ((and (listp e) (eq (first e) 'not))
                `(not-predicate ,(transform (second e))))
               ;; Already a predicate
               (t e))))
    (transform expr)))

;;; ============================================================
;;; Boundary Condition Types
;;; ============================================================

(defclass boundary-condition ()
  ((name :initarg :name
         :accessor boundary-name
         :documentation "Region/boundary name")
   (attribute :initarg :attribute
              :accessor boundary-attribute
              :initform 1
              :documentation "Mesh boundary attribute number")
   (type :initarg :type
         :accessor boundary-type
         :documentation "BC type (:velocity, :pressure, :temperature)")
   (predicate :initarg :predicate
              :accessor boundary-predicate
              :documentation "Predicate for selecting boundary elements")
   (function :initarg :function
             :accessor boundary-function
             :documentation "Function (x y [z]) -> value"))
  (:documentation "Boundary condition specification"))

(defmethod print-object ((obj boundary-condition) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~A (~A, attr=~D)"
            (boundary-name obj)
            (boundary-type obj)
            (boundary-attribute obj))))

(defun make-boundary-condition (name predicate type function)
  "Create a boundary condition."
  (make-instance 'boundary-condition
                 :name name
                 :predicate predicate
                 :type type
                 :function function))

;;; ============================================================
;;; Condition Type Constructors
;;; ============================================================

(defun velocity (&key u-function v-function w-function function value)
  "Create a velocity boundary condition specification.

   Keyword arguments:
     :u-function - Function for u-velocity component
     :v-function - Function for v-velocity component
     :w-function - Function for w-velocity component (3D)
     :function   - Single function for scalar velocity
     :value      - Constant velocity value"
  (list :type :velocity
        :u-function u-function
        :v-function v-function
        :w-function w-function
        :function (or function
                      (when value (constantly value)))))

(defun pressure (function-or-value)
  "Create a pressure boundary condition specification."
  (list :type :pressure
        :function (if (functionp function-or-value)
                      function-or-value
                      (constantly function-or-value))))

(defun temperature (function-or-value)
  "Create a temperature boundary condition specification."
  (list :type :temperature
        :function (if (functionp function-or-value)
                      function-or-value
                      (constantly function-or-value))))

;;; ============================================================
;;; Common Boundary Conditions
;;; ============================================================

(defun no-slip ()
  "Create a no-slip (zero velocity) boundary condition."
  (velocity (constantly 0.0)))

(defun inlet (velocity-profile)
  "Create an inlet boundary condition with given velocity profile."
  (velocity velocity-profile))

(defun outlet (&optional (pressure-value 0.0))
  "Create an outlet boundary condition with given pressure."
  (pressure pressure-value))

;;; ============================================================
;;; Region Macro (for SDL)
;;; ============================================================

(defmacro region (name predicate condition)
  "Define a boundary condition region.

   NAME - String name for the region
   PREDICATE - A where expression selecting boundary elements
   CONDITION - A condition specification (velocity, pressure, etc.)

   Example:
   (region \"lid\" (where (= y 1.0)) (velocity #'lid-velocity))"
  `(make-boundary-condition
    ,name
    ,predicate
    (getf ,condition :type)
    (getf ,condition :function)))
