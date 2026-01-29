;;;; demo/cavity.lisp - Lid-Driven Cavity Flow
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; This file demonstrates the SDL (Simulation Definition Language)
;;;; for defining a classic lid-driven cavity flow problem.
;;;;
;;;; Run with:
;;;;   MfemRun -f demo/cavity.lisp

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
  ;;; 10x10 quad mesh
  ;;; --------------------------------------------------------
  (mesh
    (generate *domain*
              :type :quad
              :divisions '(10 10)))

  ;;; --------------------------------------------------------
  ;;; Boundary conditions
  ;;; --------------------------------------------------------
  (boundaries
    ;; Lid moves with unit velocity
    (defun lid-velocity (x y)
      (declare (ignore x y))
      1.0)

    ;; Walls are stationary (no-slip)
    (defun wall-velocity (x y)
      (declare (ignore x y))
      0.0)

    ;; Apply BCs to each boundary
    (region "lid"    (where (= y 1.0)) (velocity #'lid-velocity))
    (region "bottom" (where (= y 0.0)) (velocity #'wall-velocity))
    (region "left"   (where (= x 0.0)) (velocity #'wall-velocity))
    (region "right"  (where (= x 1.0)) (velocity #'wall-velocity)))

  ;;; --------------------------------------------------------
  ;;; Physics: incompressible Navier-Stokes
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
