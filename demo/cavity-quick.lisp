;;;; demo/cavity-quick.lisp - Quick Test Cavity
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Coarse mesh and short runtime for quick testing.
;;;; Not suitable for validation - use cavity.lisp for that.
;;;;
;;;; Run with: StreamVorti demo/cavity-quick.lisp

(in-package :sdl)

(simulation "cavity-quick" :dim 2

  ;; Coarse mesh for speed
  (domain (box (0 0) (1 1)) :mesh :n (20 20))

  (boundaries
    (lid    (= y 1))
    (bottom (= y 0))
    (left   (= x 0))
    (right  (= x 1)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    (bc lid    :velocity (1 0))
    (bc bottom :no-slip)
    (bc left   :no-slip)
    (bc right  :no-slip))

  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  ;; Short runtime
  (temporal :explicit-euler
    :dt 0.005
    :end 2.0)

  (output :vtk
    :directory "results/quick/"
    :every 0.5
    :fields (vorticity streamfunction velocity)))
