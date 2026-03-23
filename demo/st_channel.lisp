;;;; demo/st_channel.lisp - Space-Time Channel Flow (Poiseuille)
;;;;
;;;; StreamVorti SDL Example: Space-Time Navier-Stokes
;;;;
;;;; Plane Poiseuille flow in a channel with parabolic inflow.
;;;; Exact solution: u(y) = 4 * U_max * y * (1-y), v = 0, p = -8*mu*x
;;;;
;;;; This is a useful validation case since the exact solution is known.
;;;;
;;;; Run with:
;;;;   ./SpaceTimeNS -f demo/st_channel.lisp -lp lisp -pv

(in-package :sdl)

(defparameter *channel-length* 5.0)
(defparameter *channel-height* 1.0)
(defparameter *U-max* 1.0)

(simulation "channel-flow-st" :dim 2

  (domain (box (0 0) (*channel-length* *channel-height*))
    :mesh :n (50 10)
    :element-type :quad)

  (boundaries
    (inlet  (= x 0))
    (outlet (= x *channel-length*))
    (top    (= y *channel-height*))
    (bottom (= y 0)))

  (physics :navier-stokes
    :formulation :velocity-pressure
    :method :space-time
    :Re 100

    (bc inlet  :inflow (lambda (x y t)
                         (values (* 4 *U-max* y (- *channel-height* y)) 0)))
    (bc outlet :outflow)
    (bc top    :no-slip)
    (bc bottom :no-slip))

  (spatial :fem
    :element-type :tensor-product
    :velocity-order 2
    :pressure-order 1)

  (space-time
    :formulation :st-supg-pspg
    :temporal-scheme :dg
    :dt-slab 0.5
    :stabilization (:supg :pspg :lsic))

  (nonlinear-solver
    :type :newton
    :max-iterations 15
    :rel-tolerance 1e-10)

  (linear-solver
    :type :gmres
    :preconditioner :ilu)

  (output :vtk
    :directory "results/st_channel/"
    :every 1
    :fields (velocity pressure)))
