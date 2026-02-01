;;;; demo/poiseuille.lisp - Poiseuille Channel Flow
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Pressure-driven flow between parallel plates.
;;;; Parabolic velocity profile: u(y) = 4*U_max*y*(H-y)/H²
;;;;
;;;; Run with: StreamVorti demo/poiseuille.lisp

(in-package :sdl)

(simulation "poiseuille-flow" :dim 2

  ;; Channel: length 4, height 1
  (domain (box (0 0) (4 1)) :mesh :n (80 20))

  (boundaries
    (inlet  (= x 0))
    (outlet (= x 4))
    (top    (= y 1))
    (bottom (= y 0)))

  (physics :navier-stokes
    :formulation :velocity-pressure
    :Re 100

    ;; Parabolic inlet profile
    (bc inlet :velocity
      :u (fn (x y) (* 4 y (- 1 y)))
      :v 0)
    (bc outlet :pressure 0)
    (bc top    :no-slip)
    (bc bottom :no-slip))

  (method :dcpse
    :neighbors 25
    :support-radius 5.0)

  (time :method :explicit-euler
    :dt 0.001
    :end 5.0)

  (output :vtk
    :directory "results/poiseuille/"
    :every 0.5
    :fields (velocity pressure)))
