;;;; demo/cylinder.lisp - Channel Flow (Cylinder Placeholder)
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Channel flow in a rectangular domain as a placeholder for
;;;; flow around a cylinder. The cylinder obstacle requires CSG
;;;; geometry (SDL::DIFFERENCE) and an unstructured mesh generator
;;;; (e.g. Gmsh) which are not yet supported by the C++ backend.
;;;;
;;;; Run with: StreamVorti demo/cylinder.lisp

(in-package :sdl)

(defparameter *channel-length* 5.0)
(defparameter *channel-height* 1.0)

(simulation "cylinder-flow" :dim 2

  ;; Rectangular channel: placeholder until CSG + Gmsh support is added
  (domain (box (0 0) (*channel-length* *channel-height*))
    :mesh :n (50 10))

  (boundaries
    (inlet   (= x 0))
    (outlet  (= x *channel-length*))
    (top     (= y *channel-height*))
    (bottom  (= y 0)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    (bc inlet   :velocity (1 0))
    (bc outlet  :no-slip)
    (bc top     :no-slip)
    (bc bottom  :no-slip))

  (spatial :dcpse
    :neighbors 25)

  (temporal :explicit-euler
    :dt 0.001
    :end 5.0)

  (output :vtk
    :every 1.0))
