;;;; demo/channel.lisp - Channel Flow with Parabolic Inlet
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; This file demonstrates the SDL (Simulation Definition Language)
;;;; for defining a channel flow problem with a parabolic inlet profile.
;;;;
;;;; Run with:
;;;;   MfemRun -f demo/channel.lisp

(in-package :sdl)

;;; ============================================================
;;; Simulation Definition
;;; ============================================================

(simulation
  :name "channel-flow-2d"
  :version 1
  :dimension 2

  ;;; --------------------------------------------------------
  ;;; Domain: 4:1 aspect ratio channel
  ;;; --------------------------------------------------------
  (geometry
    (defparameter *channel* (rectangle '(0 0) '(4 1))))

  ;;; --------------------------------------------------------
  ;;; Mesh: structured quadrilateral, finer in x direction
  ;;; --------------------------------------------------------
  (mesh
    (generate *channel*
              :type :quad
              :divisions '(40 10)))

  ;;; --------------------------------------------------------
  ;;; Boundary conditions
  ;;; --------------------------------------------------------
  (boundaries
    ;; Parabolic velocity profile at inlet: u(y) = 4*y*(1-y)
    ;; Maximum velocity at y=0.5, zero at walls
    (defun parabolic-inlet (x y)
      (declare (ignore x))
      (* 4.0 y (- 1.0 y)))

    ;; No-slip walls
    (defun no-slip (x y)
      (declare (ignore x y))
      0.0)

    ;; Zero pressure at outlet
    (defun outlet-pressure (x y)
      (declare (ignore x y))
      0.0)

    ;; Apply boundary conditions
    (region "inlet"  (where (= x 0.0)) (velocity #'parabolic-inlet))
    (region "outlet" (where (= x 4.0)) (pressure #'outlet-pressure))
    (region "top"    (where (= y 1.0)) (velocity #'no-slip))
    (region "bottom" (where (= y 0.0)) (velocity #'no-slip)))

  ;;; --------------------------------------------------------
  ;;; Physics
  ;;; --------------------------------------------------------
  (physics
    :type :incompressible-navier-stokes
    :formulation :stream-vorticity
    :reynolds 200.0
    :density 1.0
    :viscosity 0.005)

  ;;; --------------------------------------------------------
  ;;; Discretization
  ;;; --------------------------------------------------------
  (discretization
    :method :dcpse
    :num-neighbors 30
    :cutoff-radius 30.0
    :support-radius 5.0)

  ;;; --------------------------------------------------------
  ;;; Solver
  ;;; --------------------------------------------------------
  (solver
    :timestepping :explicit-euler
    :dt 0.0005
    :end-time 20.0
    :tolerance 1.0e-6)

  ;;; --------------------------------------------------------
  ;;; Output
  ;;; --------------------------------------------------------
  (output
    :format :vtk
    :interval 0.5
    :directory "results/channel/"
    :fields '(vorticity streamfunction velocity pressure)))

;;; ============================================================
;;; Print simulation summary
;;; ============================================================

(format t "~%Loaded: ~A~%" (simulation-name *current-simulation*))
(print-simulation-summary *current-simulation*)
