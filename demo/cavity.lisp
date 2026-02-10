;;;; demo/cavity.lisp - Lid-Driven Cavity Flow
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Classic lid-driven cavity benchmark for Re=100.
;;;; Reference: Ghia, Ghia & Shin (1982)
;;;;
;;;; Run with: StreamVorti demo/cavity.lisp

(in-package :sdl)

(simulation "lid-driven-cavity" :dim 2

  ;; Domain: structured mesh on unit square
  (domain (box (0 0) (1 1)) :mesh :n (40 40))

  ;; Named parts of ∂Ω
  (boundaries
    (lid    (= y 1))
    (bottom (= y 0))
    (left   (= x 0))
    (right  (= x 1)))

  ;; Navier-Stokes with stream-vorticity formulation
  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    (bc lid    :velocity (1 0))
    (bc bottom :no-slip)
    (bc left   :no-slip)
    (bc right  :no-slip))

  ;; DCPSE meshless discretization
  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  ;; Time integration
  (temporal :explicit-euler
    :dt 0.001
    :end 10.0
    :tolerance 1e-6)

  ;; VTK output for ParaView
  (output :vtk
    :directory "results/cavity/"
    :every 0.1
    :fields (vorticity streamfunction velocity)))
