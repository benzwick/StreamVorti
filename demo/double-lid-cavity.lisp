;;;; demo/double-lid-cavity.lisp - Double Lid-Driven Cavity
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Both top and bottom walls move in opposite directions,
;;;; creating symmetric counter-rotating vortices.
;;;;
;;;; Run with: StreamVorti demo/double-lid-cavity.lisp

(in-package :sdl)

(simulation "double-lid-cavity" :dim 2

  (domain (box (0 0) (1 1)) :mesh :n (40 40))

  (boundaries
    (top    (= y 1))
    (bottom (= y 0))
    (left   (= x 0))
    (right  (= x 1)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    ;; Top lid moves right, bottom lid moves left
    (bc top    :velocity ( 1 0))
    (bc bottom :velocity (-1 0))
    (bc left   :no-slip)
    (bc right  :no-slip))

  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  (temporal :explicit-euler
    :dt 0.001
    :end 15.0)

  (output :vtk
    :directory "results/double-lid/"
    :every 0.5
    :fields (vorticity streamfunction velocity)))
