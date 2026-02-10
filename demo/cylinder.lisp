;;;; demo/cylinder.lisp - Flow Around Cylinder
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Unsteady flow around circular cylinder.
;;;; Vortex shedding occurs for Re > 47.
;;;;
;;;; Run with: StreamVorti demo/cylinder.lisp

(in-package :sdl)

;; Domain parameters
(defparameter *domain-width* 20.0)
(defparameter *domain-height* 8.0)
(defparameter *cylinder-x* 5.0)
(defparameter *cylinder-y* 4.0)
(defparameter *cylinder-radius* 0.5)

(simulation "cylinder-flow" :dim 2

  ;; Rectangular domain with circular hole
  (domain (difference
            (box (0 0) (*domain-width* *domain-height*))
            (ball (*cylinder-x* *cylinder-y*) *cylinder-radius*))
    :mesh :generator :gmsh :h 0.2)

  (boundaries
    (inlet    (= x 0))
    (outlet   (= x *domain-width*))
    (top      (= y *domain-height*))
    (bottom   (= y 0))
    (cylinder (on-surface "cylinder")))

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
