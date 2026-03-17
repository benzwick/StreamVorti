;;;; demo/von-karman.lisp - Von Karman Vortex Street
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Flow past a circular cylinder in a channel, producing periodic
;;;; vortex shedding (von Karman vortex street) at Re ~ 100.
;;;;
;;;; Reference parameters from:
;;;;   von_karman/von_karman/cylinder/src/ns_solver.cpp
;;;;
;;;; Run with: StreamVorti -f demo/von-karman.lisp -lp lisp -pv

(in-package :sdl)

;; Geometry parameters (from von Karman reference code)
(defparameter *channel-length* 2.2d0)
(defparameter *channel-height* 0.41d0)
(defparameter *cylinder-x* 0.2d0)
(defparameter *cylinder-y* 0.2d0)
(defparameter *cylinder-r* 0.05d0)

;; Flow parameters
;; U_mean = 0.2 m/s, D = 0.1 m → Re = U_mean * D / ν = 100
(defparameter *u-mean* 0.2d0)
(defparameter *kinematic-viscosity* 0.001d0)

(simulation "von-karman-vortex-street" :dim 2

  ;; Channel with circular cylinder cut out
  (domain (difference
            (box (0 0) (*channel-length* *channel-height*))
            (ball (*cylinder-x* *cylinder-y*) *cylinder-r*))
    :mesh :generator :gmsh :h 0.01)

  (boundaries
    (inlet    (= x 0))
    (outlet   (= x *channel-length*))
    (top      (= y *channel-height*))
    (bottom   (= y 0))
    (cylinder (on-surface "cylinder")))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 100

    ;; Parabolic inlet profile: u(y) = 4 * U_mean * y * (H - y) / H²
    (bc inlet :velocity
      :u (fn (x y)
           (let ((h *channel-height*))
             (* 4 *u-mean* y (- h y) (/ 1 (* h h)))))
      :v 0)
    (bc outlet   :outflow)
    (bc top      :no-slip)
    (bc bottom   :no-slip)
    (bc cylinder :no-slip))

  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  (temporal :explicit-euler
    :dt 0.001
    :end 8.0)

  (output :vtk
    :directory "results/von-karman/"
    :every 0.1
    :fields (vorticity streamfunction velocity)))
