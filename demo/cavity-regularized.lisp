;;;; demo/cavity-regularized.lisp - Regularized Lid-Driven Cavity
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Smooth lid velocity profile: u(x) = 16*x²*(1-x)²
;;;; Zero at corners, maximum 1.0 at x=0.5.
;;;; Avoids corner singularities for better convergence.
;;;;
;;;; Run with: StreamVorti demo/cavity-regularized.lisp

(in-package :sdl)

(simulation "regularized-cavity" :dim 2

  (domain (box (0 0) (1 1)) :mesh :n (40 40))

  (boundaries
    (lid    (= y 1))
    (bottom (= y 0))
    (left   (= x 0))
    (right  (= x 1)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    ;; Regularized lid: smooth velocity profile
    (bc lid :velocity
      :u (fn (x y) (* 16 x x (- 1 x) (- 1 x)))
      :v 0)
    (bc bottom :no-slip)
    (bc left   :no-slip)
    (bc right  :no-slip))

  (method :dcpse
    :neighbors 25
    :support-radius 5.0)

  (time :method :explicit-euler
    :dt 0.001
    :end 10.0)

  (output :vtk
    :directory "results/regularized/"
    :every 0.1
    :fields (vorticity streamfunction velocity)))
