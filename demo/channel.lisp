;;;; demo/channel.lisp - Channel Flow with Parabolic Inlet
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Channel flow with parabolic inlet velocity profile using
;;;; stream function - vorticity formulation.
;;;;
;;;; Run with: StreamVorti -f demo/channel.lisp -lp lisp -pv

(in-package :sdl)

(simulation "channel-flow-2d" :dim 2

  ;; 4:1 aspect ratio channel
  (domain (box (0 0) (4 1)) :mesh :n (80 20))

  (boundaries
    (inlet  (= x 0))
    (outlet (= x 4))
    (top    (= y 1))
    (bottom (= y 0)))

  (physics :navier-stokes
    :formulation :stream-vorticity
    :Re 200

    ;; Parabolic inlet profile: u(y) = 4*y*(1-y)
    (bc inlet :velocity
      :u (fn (x y) (* 4 y (- 1 y)))
      :v 0)
    (bc outlet :pressure 0)
    (bc top    :no-slip)
    (bc bottom :no-slip))

  (spatial :dcpse
    :neighbors 25
    :support-radius 5.0)

  (temporal :explicit-euler
    :dt 0.001
    :end 5.0)

  (probes
    (line quarter-channel (x 1.0))
    (line mid-channel     (x 2.0))
    (line exit-channel    (x 3.0)))

  (output :vtk
    :directory "results/channel/"
    :every 1.0
    :fields (vorticity streamfunction velocity)))
