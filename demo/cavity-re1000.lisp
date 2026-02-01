;;;; demo/cavity-re1000.lisp - Lid-Driven Cavity at Re=1000
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Higher Reynolds number with prominent secondary corner vortices.
;;;; Reference: Ghia, Ghia & Shin (1982)
;;;;
;;;; Run with: StreamVorti demo/cavity-re1000.lisp

(in-package :sdl)

(simulation "lid-driven-cavity-re1000" :dim 2

  ;; Finer mesh for higher Re
  (domain (box (0 0) (1 1)) :mesh :n (60 60))

  (boundaries
    (lid    (= y 1))
    (bottom (= y 0))
    (left   (= x 0))
    (right  (= x 1)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 1000

    (bc lid    :velocity (1 0))
    (bc bottom :no-slip)
    (bc left   :no-slip)
    (bc right  :no-slip))

  (method :dcpse
    :neighbors 25
    :support-radius 5.0)

  ;; Smaller timestep for stability
  (time :method :explicit-euler
    :dt 0.0005
    :end 30.0
    :tolerance 1e-6)

  (output :vtk
    :directory "results/cavity-re1000/"
    :every 1.0
    :fields (vorticity streamfunction velocity)))
