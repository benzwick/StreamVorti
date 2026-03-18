;;;; demo/cavity-fdm.lisp - Lid-Driven Cavity Flow with Finite Differences
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Same benchmark as cavity.lisp but using finite difference (FDM)
;;;; spatial discretization instead of DCPSE.
;;;; Reference: Ghia, Ghia & Shin (1982)
;;;;
;;;; Run with: StreamVorti -f demo/cavity-fdm.lisp -pv

(in-package :sdl)

(simulation "lid-driven-cavity-fdm" :dim 2

  ;; Domain: structured mesh on unit square
  (domain (box (0 0) (1 1)) :mesh :n (40 40))

  ;; Named parts of boundary
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

  ;; Finite difference discretization (structured grid required)
  (spatial :fdm)

  ;; Time integration
  (temporal :explicit-euler
    :dt 0.001
    :end 10.0
    :tolerance 1e-6)

  ;; Line probes for Ghia validation
  (probes
    (line vertical-center   (x 0.5))
    (line horizontal-center (y 0.5)))

  ;; VTK output for ParaView
  (output :vtk
    :directory "results/cavity-fdm/"
    :every 0.1
    :fields (vorticity streamfunction velocity)))
