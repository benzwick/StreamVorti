;;;; demo/cylinder.lisp - Flow Around Cylinder
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Unsteady flow around circular cylinder.
;;;; Vortex shedding occurs for Re > 47.
;;;;
;;;; Requires: STREAMVORTI_WITH_GMSH=ON
;;;; Run with: StreamVorti -f demo/cylinder.lisp -lp lisp

(in-package :sdl)

;; Domain parameters
(defparameter *L* 20.0d0)   ; channel length
(defparameter *H* 8.0d0)    ; channel height
(defparameter *cx* 5.0d0)   ; cylinder center x
(defparameter *cy* 4.0d0)   ; cylinder center y
(defparameter *r* 0.5d0)    ; cylinder radius

(simulation "cylinder-flow" :dim 2

  (domain :gmsh
    ;; OCC geometry: rectangle with circular hole
    (occ:rectangle 0 0 0 *L* *H* :tag 1)
    (occ:disk *cx* *cy* 0 *r* *r* :tag 2)
    (occ:cut (gmsh:surface-tags '(1)) (gmsh:surface-tags '(2)))
    (occ:synchronize)

    ;; Identify boundaries (see Gmsh tutorial t16):
    ;; get-boundary returns all boundary curves of the result,
    ;; get-closest-entities selects curves nearest to a given point.
    (let ((bnd (gmsh:get-boundary (gmsh:get-entities :dim 2) :oriented nil)))
      (gmsh:add-physical-group 1
        (gmsh:tags-of (occ:get-closest-entities 0 (/ *H* 2) 0 bnd :n 1))
        :tag 1 :name "inlet")
      (gmsh:add-physical-group 1
        (gmsh:tags-of (occ:get-closest-entities *L* (/ *H* 2) 0 bnd :n 1))
        :tag 2 :name "outlet")
      (gmsh:add-physical-group 1
        (gmsh:tags-of (occ:get-closest-entities (/ *L* 2) *H* 0 bnd :n 1))
        :tag 3 :name "top")
      (gmsh:add-physical-group 1
        (gmsh:tags-of (occ:get-closest-entities (/ *L* 2) 0 0 bnd :n 1))
        :tag 4 :name "bottom")
      (gmsh:add-physical-group 1
        (gmsh:tags-of (occ:get-closest-entities *cx* *cy* 0 bnd :n 4))
        :tag 5 :name "cylinder"))
    (gmsh:add-physical-group 2
      (gmsh:tags-of (gmsh:get-entities :dim 2))
      :tag 1 :name "fluid")

    ;; Mesh refinement near cylinder and wake
    (mesh:field "Distance" :tag 1
      :number-lists '(("CurvesList" (5 6 7 8)))
      :numbers '(("Sampling" 100)))
    (mesh:field "Threshold" :tag 2
      :numbers `(("InField" 1) ("SizeMin" 0.1) ("SizeMax" 0.5)
                 ("DistMin" 0) ("DistMax" 2)))
    (mesh:field "Box" :tag 3
      :numbers `(("VIn" 0.2) ("VOut" 0.5)
                 ("XMin" ,*cx*) ("XMax" ,*L*)
                 ("YMin" ,(- *cy* (* 3 *r*)))
                 ("YMax" ,(+ *cy* (* 3 *r*)))))
    (mesh:field "Min" :tag 4
      :number-lists '(("FieldsList" (2 3)))
      :as-background t)

    (gmsh/mesh:generate :dim 2))

  ;; No (boundaries ...) needed — derived from physical group names above.

  (physics :navier-stokes
    :formulation :velocity-pressure
    :Re 100

    (bc inlet    :velocity (1 0))
    (bc outlet   :outflow)
    (bc top      :slip)
    (bc bottom   :slip)
    (bc cylinder :no-slip))

  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  (temporal :explicit-euler
    :dt 0.005
    :end 50.0)

  (output :vtk
    :directory "results/cylinder/"
    :every 0.5
    :fields (velocity pressure vorticity)))
