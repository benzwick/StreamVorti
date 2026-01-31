;;;; mesh.lisp - Mesh generation and manipulation
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :streamvorti.mesh)

;;; ============================================================
;;; Mesh Class
;;; ============================================================

(defclass mesh ()
  ((ptr :initarg :ptr
        :accessor mesh-ptr
        :documentation "Pointer to MFEM Mesh object")
   (dimension :initarg :dimension
              :accessor mesh-dimension
              :documentation "Mesh dimension (2 or 3)")
   (geometry :initarg :geometry
             :accessor mesh-geometry
             :initform nil
             :documentation "Source geometry shape (if generated)")
   (element-type :initarg :element-type
                 :accessor mesh-element-type
                 :documentation "Element type (:quad, :tri, :hex, :tet)")
   (divisions :initarg :divisions
              :accessor mesh-divisions
              :initform nil
              :documentation "Number of divisions per direction"))
  (:documentation "Wrapper for MFEM Mesh"))

(defmethod mesh-num-vertices ((m mesh))
  "Get number of vertices in mesh."
  (streamvorti.ffi:mesh-num-vertices (mesh-ptr m)))

(defmethod mesh-num-elements ((m mesh))
  "Get number of elements in mesh."
  (streamvorti.ffi:mesh-num-elements (mesh-ptr m)))

(defmethod print-object ((obj mesh) stream)
  (print-unreadable-object (obj stream :type t)
    (format stream "~DD, ~A vertices, ~A elements"
            (mesh-dimension obj)
            (mesh-num-vertices obj)
            (mesh-num-elements obj))))

;;; ============================================================
;;; Mesh Generation
;;; ============================================================

(defun element-type-code (type-keyword dimension)
  "Convert element type keyword to MFEM element type code."
  (cond
    ((and (= dimension 2) (member type-keyword '(:quad :quadrilateral)))
     3)
    ((and (= dimension 2) (member type-keyword '(:tri :triangle)))
     2)
    ((and (= dimension 3) (member type-keyword '(:hex :hexahedron)))
     5)
    ((and (= dimension 3) (member type-keyword '(:tet :tetrahedron)))
     4)
    (t (error "Unknown element type ~A for dimension ~D" type-keyword dimension))))

(defun generate-mesh (geometry &key (type :quad) (divisions '(10 10))
                                    (size-x 1.0) (size-y 1.0) (size-z 1.0))
  "Generate a mesh from geometry specification.

   GEOMETRY - A shape object (rectangle, box, etc.) or NIL for unit domain
   TYPE - Element type: :quad, :tri (2D) or :hex, :tet (3D)
   DIVISIONS - List of (nx ny) or (nx ny nz) divisions
   SIZE-X/Y/Z - Domain size (used if geometry is NIL)"
  (let ((dim (cond
               ((typep geometry 'streamvorti.geometry:rectangle) 2)
               ((typep geometry 'streamvorti.geometry:box) 3)
               ((= (length divisions) 2) 2)
               ((= (length divisions) 3) 3)
               (t 2))))

    ;; Get sizes from geometry if available
    (when geometry
      (let ((bounds (streamvorti.geometry:shape-bounds geometry)))
        (let ((min-pt (first bounds))
              (max-pt (second bounds)))
          (setf size-x (- (first max-pt) (first min-pt)))
          (setf size-y (- (second max-pt) (second min-pt)))
          (when (= dim 3)
            (setf size-z (- (third max-pt) (third min-pt)))))))

    (let* ((nx (or (first divisions) 10))
           (ny (or (second divisions) 10))
           (nz (or (third divisions) 10))
           (elem-code (element-type-code type dim))
           (ptr (if (= dim 2)
                    (streamvorti.ffi:make-cartesian-mesh-2d nx ny elem-code
                                                            size-x size-y)
                    (streamvorti.ffi:make-cartesian-mesh-3d nx ny nz elem-code
                                                            size-x size-y size-z))))
      (make-instance 'mesh
                     :ptr ptr
                     :dimension dim
                     :geometry geometry
                     :element-type type
                     :divisions divisions))))

(defun load-mesh-file (path &key (format :auto))
  "Load a mesh from file.

   PATH - File path
   FORMAT - File format (:mfem, :gmsh, :vtk, :netgen, or :auto)"
  (declare (ignore format))
  (let ((ptr (streamvorti.ffi:load-mesh path 0)))
    (unless ptr
      (error "Failed to load mesh from ~A" path))
    (let ((dim (streamvorti.ffi:mesh-dimension ptr)))
      (make-instance 'mesh
                     :ptr ptr
                     :dimension dim
                     :element-type (if (= dim 2) :quad :hex)))))

(defun save-mesh-file (mesh path &key (format :mfem))
  "Save mesh to file.

   MESH - Mesh object
   PATH - Output file path
   FORMAT - Output format (:mfem or :vtk)"
  (let ((format-code (case format
                       (:vtk 2)
                       (otherwise 1))))
    (streamvorti.ffi:save-mesh (mesh-ptr mesh) path format-code)))

(defun refine-mesh (mesh &key (times 1))
  "Uniformly refine mesh.

   MESH - Mesh object
   TIMES - Number of refinement iterations"
  (streamvorti.ffi:mesh-refine (mesh-ptr mesh) times)
  mesh)

;;; ============================================================
;;; Mesh Specification (for SDL)
;;; ============================================================

(defclass mesh-spec ()
  ((type :initarg :type
         :accessor mesh-spec-type
         :documentation "Type: :generated or :loaded")
   (path :initarg :path
         :accessor mesh-spec-path
         :initform nil
         :documentation "File path (for loaded meshes)")
   (geometry :initarg :geometry
             :accessor mesh-spec-geometry
             :initform nil
             :documentation "Geometry shape (for generated meshes)")
   (element-type :initarg :element-type
                 :accessor mesh-spec-element-type
                 :initform :quad
                 :documentation "Element type keyword")
   (divisions :initarg :divisions
              :accessor mesh-spec-divisions
              :initform '(10 10)
              :documentation "Mesh divisions")
   (sizes :initarg :sizes
          :accessor mesh-spec-sizes
          :initform '(1.0 1.0)
          :documentation "Domain sizes"))
  (:documentation "Specification for mesh generation/loading"))

(defun make-mesh-spec-generate (geometry &key (type :quad) (divisions '(10 10)))
  "Create a mesh specification for generated mesh."
  (make-instance 'mesh-spec
                 :type :generated
                 :geometry geometry
                 :element-type type
                 :divisions divisions))

(defun make-mesh-spec-load (path &key (format :auto))
  "Create a mesh specification for loading mesh from file."
  (declare (ignore format))
  (make-instance 'mesh-spec
                 :type :loaded
                 :path path))

(defun realize-mesh-spec (spec)
  "Convert mesh specification to actual mesh."
  (ecase (mesh-spec-type spec)
    (:generated
     (generate-mesh (mesh-spec-geometry spec)
                    :type (mesh-spec-element-type spec)
                    :divisions (mesh-spec-divisions spec)))
    (:loaded
     (load-mesh-file (mesh-spec-path spec)))))
