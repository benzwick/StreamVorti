;;;; geometry.lisp - Geometry primitives for SDL
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :streamvorti.geometry)

;;; ============================================================
;;; Base Shape Class
;;; ============================================================

(defclass shape ()
  ((dimension :initarg :dimension
              :accessor shape-dimension
              :documentation "Spatial dimension of the shape (2 or 3)"))
  (:documentation "Base class for all geometric shapes"))

(defgeneric shape-bounds (shape)
  (:documentation "Return bounding box as (min-point max-point)"))

(defgeneric shape-contains-p (shape point)
  (:documentation "Test if shape contains point"))

;;; ============================================================
;;; 1D Shapes
;;; ============================================================

(defclass interval (shape)
  ((min-corner :initarg :min-corner
               :accessor interval-min-corner
               :documentation "Minimum point (x)")
   (max-corner :initarg :max-corner
               :accessor interval-max-corner
               :documentation "Maximum point (x)"))
  (:default-initargs :dimension 1)
  (:documentation "1D interval region"))

(defmethod shape-bounds ((shape interval))
  (list (interval-min-corner shape)
        (interval-max-corner shape)))

(defmethod shape-contains-p ((shape interval) point)
  (let ((x (if (listp point) (first point) point))
        (min-x (first (interval-min-corner shape)))
        (max-x (first (interval-max-corner shape))))
    (and (>= x min-x) (<= x max-x))))

(defmethod print-object ((obj interval) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~A to ~A"
            (interval-min-corner obj)
            (interval-max-corner obj))))

;;; ============================================================
;;; 2D Shapes
;;; ============================================================

(defclass rectangle (shape)
  ((min-corner :initarg :min-corner
               :accessor rectangle-min-corner
               :documentation "Lower-left corner (x y)")
   (max-corner :initarg :max-corner
               :accessor rectangle-max-corner
               :documentation "Upper-right corner (x y)"))
  (:default-initargs :dimension 2)
  (:documentation "2D rectangular region"))

(defun make-rectangle (min-corner max-corner)
  "Create a rectangle from min and max corners.
   Each corner is a list (x y)."
  (make-instance 'rectangle
                 :min-corner (coerce min-corner 'list)
                 :max-corner (coerce max-corner 'list)))

(defun rectangle (min-corner max-corner)
  "Alias for make-rectangle."
  (make-rectangle min-corner max-corner))

(defmethod shape-bounds ((shape rectangle))
  (list (rectangle-min-corner shape)
        (rectangle-max-corner shape)))

(defmethod shape-contains-p ((shape rectangle) point)
  (let ((min (rectangle-min-corner shape))
        (max (rectangle-max-corner shape))
        (x (first point))
        (y (second point)))
    (and (>= x (first min)) (<= x (first max))
         (>= y (second min)) (<= y (second max)))))

(defmethod print-object ((obj rectangle) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~A to ~A"
            (rectangle-min-corner obj)
            (rectangle-max-corner obj))))

;;; --------------------------------------------------------

(defclass circle (shape)
  ((center :initarg :center
           :accessor circle-center
           :documentation "Center point (x y)")
   (radius :initarg :radius
           :accessor circle-radius
           :documentation "Circle radius"))
  (:default-initargs :dimension 2)
  (:documentation "2D circular region"))

(defun make-circle (center radius)
  "Create a circle from center point and radius."
  (make-instance 'circle
                 :center (coerce center 'list)
                 :radius (coerce radius 'double-float)))

(defun circle (center radius)
  "Alias for make-circle."
  (make-circle center radius))

(defmethod shape-bounds ((shape circle))
  (let ((cx (first (circle-center shape)))
        (cy (second (circle-center shape)))
        (r (circle-radius shape)))
    (list (list (- cx r) (- cy r))
          (list (+ cx r) (+ cy r)))))

(defmethod shape-contains-p ((shape circle) point)
  (let ((cx (first (circle-center shape)))
        (cy (second (circle-center shape)))
        (r (circle-radius shape))
        (x (first point))
        (y (second point)))
    (<= (+ (expt (- x cx) 2)
           (expt (- y cy) 2))
        (expt r 2))))

;;; --------------------------------------------------------

(defclass polygon (shape)
  ((vertices :initarg :vertices
             :accessor polygon-vertices
             :documentation "List of vertices ((x1 y1) (x2 y2) ...)"))
  (:default-initargs :dimension 2)
  (:documentation "2D polygonal region"))

(defun make-polygon (vertices)
  "Create a polygon from a list of vertices."
  (make-instance 'polygon
                 :vertices (mapcar (lambda (v) (coerce v 'list)) vertices)))

(defun polygon (&rest vertices)
  "Create a polygon from vertices."
  (make-polygon vertices))

(defmethod shape-bounds ((shape polygon))
  (let ((vertices (polygon-vertices shape)))
    (let ((xs (mapcar #'first vertices))
          (ys (mapcar #'second vertices)))
      (list (list (apply #'min xs) (apply #'min ys))
            (list (apply #'max xs) (apply #'max ys))))))

;;; ============================================================
;;; 3D Shapes
;;; ============================================================

(defclass box (shape)
  ((min-corner :initarg :min-corner
               :accessor box-min-corner
               :documentation "Minimum corner (x y z)")
   (max-corner :initarg :max-corner
               :accessor box-max-corner
               :documentation "Maximum corner (x y z)"))
  (:default-initargs :dimension 3)
  (:documentation "3D box/cuboid region"))

(defun make-box (min-corner max-corner)
  "Create a box from min and max corners.
   Each corner is a list (x y z)."
  (make-instance 'box
                 :min-corner (coerce min-corner 'list)
                 :max-corner (coerce max-corner 'list)))

(defun box (min-corner max-corner)
  "Alias for make-box."
  (make-box min-corner max-corner))

(defmethod shape-bounds ((shape box))
  (list (box-min-corner shape)
        (box-max-corner shape)))

(defmethod shape-contains-p ((shape box) point)
  (let ((min (box-min-corner shape))
        (max (box-max-corner shape))
        (x (first point))
        (y (second point))
        (z (third point)))
    (and (>= x (first min)) (<= x (first max))
         (>= y (second min)) (<= y (second max))
         (>= z (third min)) (<= z (third max)))))

;;; --------------------------------------------------------

(defclass sphere (shape)
  ((center :initarg :center
           :accessor sphere-center
           :documentation "Center point (x y z)")
   (radius :initarg :radius
           :accessor sphere-radius
           :documentation "Sphere radius"))
  (:default-initargs :dimension 3)
  (:documentation "3D spherical region"))

(defun make-sphere (center radius)
  "Create a sphere from center point and radius."
  (make-instance 'sphere
                 :center (coerce center 'list)
                 :radius (coerce radius 'double-float)))

(defun sphere (center radius)
  "Alias for make-sphere."
  (make-sphere center radius))

(defmethod shape-bounds ((shape sphere))
  (let ((cx (first (sphere-center shape)))
        (cy (second (sphere-center shape)))
        (cz (third (sphere-center shape)))
        (r (sphere-radius shape)))
    (list (list (- cx r) (- cy r) (- cz r))
          (list (+ cx r) (+ cy r) (+ cz r)))))

(defmethod shape-contains-p ((shape sphere) point)
  (let ((center (sphere-center shape))
        (r (sphere-radius shape)))
    (<= (+ (expt (- (first point) (first center)) 2)
           (expt (- (second point) (second center)) 2)
           (expt (- (third point) (third center)) 2))
        (expt r 2))))

;;; --------------------------------------------------------

(defclass cylinder (shape)
  ((base-center :initarg :base-center
                :accessor cylinder-base-center
                :documentation "Center of base (x y z)")
   (top-center :initarg :top-center
               :accessor cylinder-top-center
               :documentation "Center of top (x y z)")
   (radius :initarg :radius
           :accessor cylinder-radius
           :documentation "Cylinder radius"))
  (:default-initargs :dimension 3)
  (:documentation "3D cylindrical region"))

(defun make-cylinder (base-center top-center radius)
  "Create a cylinder from base center, top center, and radius."
  (make-instance 'cylinder
                 :base-center (coerce base-center 'list)
                 :top-center (coerce top-center 'list)
                 :radius (coerce radius 'double-float)))

(defun cylinder (base-center top-center radius)
  "Alias for make-cylinder."
  (make-cylinder base-center top-center radius))

;;; ============================================================
;;; CSG Operations (Future Implementation)
;;; ============================================================

(defclass csg-shape (shape)
  ((operands :initarg :operands
             :accessor csg-operands
             :documentation "List of operand shapes"))
  (:documentation "Base class for CSG operations"))

(defclass union-shape (csg-shape)
  ()
  (:documentation "CSG union of shapes"))

(defclass difference-shape (csg-shape)
  ()
  (:documentation "CSG difference of shapes"))

(defclass intersection-shape (csg-shape)
  ()
  (:documentation "CSG intersection of shapes"))

;; Constructor stubs for future implementation
(defun csg-union (&rest shapes)
  "Create a CSG union of shapes (future)."
  (make-instance 'union-shape
                 :operands shapes
                 :dimension (shape-dimension (first shapes))))

(defun csg-difference (shape &rest subtracted)
  "Create a CSG difference (future)."
  (make-instance 'difference-shape
                 :operands (cons shape subtracted)
                 :dimension (shape-dimension shape)))

(defun csg-intersection (&rest shapes)
  "Create a CSG intersection (future)."
  (make-instance 'intersection-shape
                 :operands shapes
                 :dimension (shape-dimension (first shapes))))
