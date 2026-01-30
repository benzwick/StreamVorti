;;;; demo/cavity.lisp - Lid-Driven Cavity Flow
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; This file demonstrates the SDL (Simulation Definition Language)
;;;; for defining a classic lid-driven cavity flow problem.
;;;;
;;;; Run with:
;;;;   MfemRun -f demo/cavity.lisp -pv
;;;;
;;;; For validation against Ghia et al. (1982), use a 40x40 grid.

(in-package :sdl)

;;; ============================================================
;;; Simulation Definition
;;; ============================================================

(simulation
  :name "lid-driven-cavity-2d"
  :version 1
  :dimension 2

  ;;; --------------------------------------------------------
  ;;; Domain: unit square
  ;;; --------------------------------------------------------
  (geometry
    (defparameter *domain* (rectangle '(0 0) '(1 1))))

  ;;; --------------------------------------------------------
  ;;; 40x40 quad mesh (recommended for Ghia validation)
  ;;; Use 10x10 for quick tests, 40x40 for accurate validation
  ;;; --------------------------------------------------------
  (mesh
    (generate *domain*
              :type :quad
              :divisions '(40 40)))

  ;;; --------------------------------------------------------
  ;;; Boundary conditions
  ;;; For velocity BCs, we specify separate u and v functions
  ;;; --------------------------------------------------------
  (boundaries
    ;; Lid: u=1, v=0 (moving in x-direction)
    (defun lid-u-velocity (x y)
      "Lid moves with unit velocity in x-direction"
      (declare (ignore x y))
      1.0d0)

    (defun lid-v-velocity (x y)
      "Lid has zero velocity in y-direction"
      (declare (ignore x y))
      0.0d0)

    ;; Walls: u=0, v=0 (no-slip)
    (defun wall-u-velocity (x y)
      "Walls have zero u-velocity"
      (declare (ignore x y))
      0.0d0)

    (defun wall-v-velocity (x y)
      "Walls have zero v-velocity"
      (declare (ignore x y))
      0.0d0)

    ;; Apply BCs to each boundary with separate u and v functions
    (region "lid"
      (where (= y 1.0))
      (velocity :u-function #'lid-u-velocity
                :v-function #'lid-v-velocity))

    (region "bottom"
      (where (= y 0.0))
      (velocity :u-function #'wall-u-velocity
                :v-function #'wall-v-velocity))

    (region "left"
      (where (= x 0.0))
      (velocity :u-function #'wall-u-velocity
                :v-function #'wall-v-velocity))

    (region "right"
      (where (= x 1.0))
      (velocity :u-function #'wall-u-velocity
                :v-function #'wall-v-velocity)))

  ;;; --------------------------------------------------------
  ;;; Physics: incompressible Navier-Stokes
  ;;; Re=100 for initial validation (Ghia 1982 benchmark)
  ;;; --------------------------------------------------------
  (physics
    :type :incompressible-navier-stokes
    :formulation :stream-vorticity
    :reynolds 100.0
    :density 1.0
    :viscosity 0.01)

  ;;; --------------------------------------------------------
  ;;; DCPSE discretization
  ;;; --------------------------------------------------------
  (discretization
    :method :dcpse
    :num-neighbors 25
    :cutoff-radius 30.0
    :support-radius 5.0)

  ;;; --------------------------------------------------------
  ;;; Time integration
  ;;; --------------------------------------------------------
  (solver
    :timestepping :explicit-euler
    :dt 0.001
    :end-time 10.0
    :tolerance 1.0e-6
    :max-iterations 10000)

  ;;; --------------------------------------------------------
  ;;; Output configuration
  ;;; --------------------------------------------------------
  (output
    :format :vtk
    :interval 0.1
    :directory "results/cavity/"
    :fields '(vorticity streamfunction velocity)))

;;; ============================================================
;;; Print simulation summary
;;; ============================================================

(format t "~%Loaded: ~A~%" (simulation-name *current-simulation*))
(print-simulation-summary *current-simulation*)
