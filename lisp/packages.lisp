;;;; packages.lisp - Package definitions for StreamVorti SDL
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :cl-user)

;;; ============================================================
;;; Foreign Function Interface Package
;;; ============================================================

(defpackage :streamvorti.ffi
  (:use :cl)
  (:documentation "Foreign function interface to StreamVorti C++ code")
  (:export
   ;; Mesh functions
   #:make-cartesian-mesh-2d
   #:make-cartesian-mesh-3d
   #:load-mesh
   #:save-mesh
   #:free-mesh
   #:mesh-num-vertices
   #:mesh-num-elements
   #:mesh-num-boundary-elements
   #:mesh-dimension
   #:mesh-space-dimension
   #:mesh-get-vertices
   #:mesh-set-boundary-attribute
   #:mesh-refine
   ;; GridFunction
   #:make-grid-function
   #:free-grid-function
   #:grid-function-size
   #:grid-function-get-values
   #:grid-function-set-values
   ;; DCPSE functions
   #:make-dcpse-2d
   #:make-dcpse-3d
   #:free-dcpse
   #:dcpse-update
   #:dcpse-apply-dx
   #:dcpse-apply-dy
   #:dcpse-apply-dz
   #:dcpse-apply-dxx
   #:dcpse-apply-dyy
   #:dcpse-apply-dzz
   #:dcpse-apply-dxy
   #:dcpse-apply-laplacian
   #:dcpse-save-matrix
   #:dcpse-save-neighbors))

;;; ============================================================
;;; Geometry Package
;;; ============================================================

(defpackage :streamvorti.geometry
  (:use :cl)
  (:documentation "Geometry primitives and CSG operations")
  (:export
   ;; Shape classes
   #:shape
   #:shape-dimension
   #:shape-bounds
   #:shape-contains-p
   ;; 2D primitives
   #:rectangle
   #:circle
   #:polygon
   ;; 3D primitives
   #:box
   #:sphere
   #:cylinder
   ;; CSG operations (future)
   #:union-shape
   #:difference-shape
   #:intersection-shape
   ;; Constructors
   #:make-rectangle
   #:make-circle
   #:make-box
   #:make-sphere))

;;; ============================================================
;;; Mesh Package
;;; ============================================================

(defpackage :streamvorti.mesh
  (:use :cl :streamvorti.ffi :streamvorti.geometry)
  (:documentation "Mesh generation and manipulation")
  (:export
   #:mesh
   #:mesh-ptr
   #:mesh-dimension
   #:mesh-num-vertices
   #:mesh-num-elements
   #:generate-mesh
   #:load-mesh-file
   #:save-mesh-file
   #:refine-mesh))

;;; ============================================================
;;; Boundaries Package
;;; ============================================================

(defpackage :streamvorti.boundaries
  (:use :cl)
  (:documentation "Boundary condition handling")
  (:export
   #:boundary-condition
   #:boundary-name
   #:boundary-attribute
   #:boundary-type
   #:boundary-function
   ;; Predicates
   #:where
   #:x-equals
   #:y-equals
   #:z-equals
   #:and-predicate
   #:or-predicate
   #:not-predicate
   ;; Condition types
   #:velocity
   #:pressure
   #:temperature
   ;; Common BCs
   #:no-slip
   #:inlet
   #:outlet))

;;; ============================================================
;;; DCPSE Package
;;; ============================================================

(defpackage :streamvorti.dcpse
  (:use :cl :streamvorti.ffi :streamvorti.mesh)
  (:documentation "DC-PSE discretization")
  (:export
   #:dcpse
   #:dcpse-ptr
   #:dcpse-dimension
   #:dcpse-num-neighbors
   #:make-dcpse
   #:update-dcpse
   #:apply-derivative
   #:apply-laplacian
   #:save-derivative-matrix))

;;; ============================================================
;;; Main SDL Package
;;; ============================================================

(defpackage :sdl
  (:use :cl
        :streamvorti.ffi
        :streamvorti.geometry
        :streamvorti.mesh
        :streamvorti.boundaries
        :streamvorti.dcpse)
  (:documentation "Simulation Definition Language for StreamVorti")
  (:export
   ;; Main macros
   #:simulation
   #:geometry
   #:mesh
   #:boundaries
   #:physics
   #:discretization
   #:solver
   #:output
   ;; Geometry re-exports
   #:rectangle
   #:circle
   #:polygon
   #:box
   #:sphere
   #:cylinder
   ;; Mesh
   #:generate
   #:load-mesh
   ;; Boundaries
   #:region
   #:where
   #:velocity
   #:pressure
   #:no-slip
   ;; Simulation access
   #:*current-simulation*
   #:get-current-simulation
   #:simulation-name
   #:simulation-dimension
   #:simulation-mesh
   #:simulation-boundaries
   #:simulation-physics
   #:simulation-discretization
   #:simulation-solver
   #:simulation-output
   ;; Accessors for C++ loader
   #:get-name
   #:get-version
   #:get-dimension
   #:get-mesh
   #:get-boundaries
   #:get-physics
   #:get-discretization
   #:get-solver
   #:get-output))
