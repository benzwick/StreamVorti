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
   ;; 1D primitives
   #:interval
   #:interval-min-corner
   #:interval-max-corner
   ;; 2D primitives
   #:rectangle
   #:rectangle-min-corner
   #:rectangle-max-corner
   #:circle
   #:circle-center
   #:circle-radius
   #:polygon
   ;; 3D primitives
   #:box
   #:box-min-corner
   #:box-max-corner
   #:sphere
   #:sphere-center
   #:sphere-radius
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
   #:refine-mesh
   ;; Mesh specification (for SDL)
   #:mesh-spec
   #:mesh-spec-type
   #:mesh-spec-path
   #:mesh-spec-geometry
   #:mesh-spec-element-type
   #:mesh-spec-divisions
   #:mesh-spec-sizes
   #:make-mesh-spec-generate
   #:make-mesh-spec-load
   #:realize-mesh-spec))

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
   #:boundary-u-function
   #:boundary-v-function
   #:boundary-w-function
   #:boundary-predicate
   #:predicate-test-function
   #:predicate-axis
   #:predicate-value
   #:predicate-tolerance
   ;; Predicates
   #:where
   #:x-equals
   #:y-equals
   #:z-equals
   #:x-greater-than
   #:y-greater-than
   #:z-greater-than
   #:x-less-than
   #:y-less-than
   #:z-less-than
   #:x-greater-equal
   #:y-greater-equal
   #:z-greater-equal
   #:x-less-equal
   #:y-less-equal
   #:z-less-equal
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
   #:outlet
   ;; Region macro for SDL
   #:region))

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
;;; Main SDL Package (v2)
;;; ============================================================

(defpackage :sdl
  (:use :cl)
  (:import-from :streamvorti.geometry
   #:shape #:shape-dimension #:shape-contains-p #:shape-bounds
   #:interval #:interval-min-corner #:interval-max-corner
   #:rectangle #:rectangle-min-corner #:rectangle-max-corner
   #:circle #:circle-center #:circle-radius
   #:sphere #:sphere-center #:sphere-radius
   #:make-rectangle #:make-circle
   #:make-box #:make-sphere)
  (:documentation "Simulation Definition Language v2 for StreamVorti")
  (:export
   ;; Geometry wrappers
   #:box
   #:ball
   #:geometry-dimension
   #:geometry-min
   #:geometry-max
   #:geometry-contains-p
   #:ball-center
   #:ball-radius
   ;; Domain
   #:domain
   #:domain-type
   #:domain-dimension
   #:domain-n
   #:domain-h
   #:domain-file
   #:domain-particle-count
   ;; Predicates
   #:make-predicate
   #:predicate-matches-p
   ;; Boundary selectors
   #:make-boundary-selector
   #:selector-type
   #:selector-attribute
   ;; Boundaries
   #:make-boundaries
   #:boundary-name
   ;; Subdomains
   #:make-subdomains
   #:subdomain-name
   ;; Physics
   #:make-physics
   #:physics-type
   #:physics-Re
   #:physics-conductivity
   #:physics-formulation
   #:physics-bcs
   ;; Boundary conditions
   #:make-bc
   #:bc-type
   #:bc-boundary
   #:bc-value
   #:bc-function
   #:bc-h
   #:bc-T-inf
   ;; Spatial discretization
   #:make-spatial
   #:spatial-type
   #:spatial-neighbors
   #:spatial-support-radius
   #:spatial-order
   #:spatial-kernel
   #:spatial-h
   ;; Temporal integration
   #:make-temporal
   #:temporal-method
   #:temporal-dt
   #:temporal-end
   #:temporal-tolerance
   #:make-steady-solver
   #:make-linear-solver
   #:linear-solver-method
   #:linear-solver-preconditioner
   ;; Materials
   #:make-material
   #:material-name
   #:material-density
   #:material-youngs-modulus
   #:material-poissons-ratio
   #:material-viscosity
   #:material-conductivity
   ;; Output
   #:make-output
   #:output-format
   #:output-directory
   #:output-interval
   #:output-fields
   #:output-probes
   ;; Coupling
   #:make-coupling
   #:coupling-name
   #:coupling-type
   #:coupling-physics
   #:coupling-iterations
   #:coupling-interface
   ;; Simulation
   #:simulation
   #:*current-simulation*
   #:get-current-simulation
   #:simulation-name
   #:simulation-dim
   #:simulation-domain
   #:simulation-boundaries
   #:simulation-physics
   #:simulation-subdomains
   #:simulation-physics-list
   #:simulation-spatial
   #:simulation-temporal
   #:simulation-coupling
   ;; C++ bridge accessors (used by Loader::getProperty)
   #:get-name
   #:get-version
   #:get-dimension
   #:get-mesh
   #:get-boundaries
   #:get-physics
   #:get-discretization
   #:get-solver
   #:get-output
   #:get-type
   #:get-divisions
   #:get-size-x
   #:get-size-y
   #:get-size-z
   #:get-element-type
   #:get-path
   #:get-attribute
   #:get-predicate-axis
   #:get-predicate-value
   #:get-predicate-tolerance
   #:get-function
   #:get-u-function
   #:get-v-function
   #:get-w-function
   #:get-formulation
   #:get-reynolds
   #:get-density
   #:get-viscosity
   #:get-num-neighbors
   #:get-cutoff-radius
   #:get-support-radius
   #:get-timestepping
   #:get-dt
   #:get-end-time
   #:get-tolerance
   #:get-max-iterations
   #:get-format
   #:get-interval
   #:get-directory
   #:get-fields))
