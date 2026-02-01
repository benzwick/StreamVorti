;;;; demo/backward-step.lisp - Backward Facing Step
;;;;
;;;; StreamVorti SDL Example
;;;;
;;;; Classic benchmark for separated flow and reattachment.
;;;; Reference: Armaly et al. (1983)
;;;;
;;;; Run with: StreamVorti demo/backward-step.lisp

(in-package :sdl)

;; Geometry parameters (used in multiple places)
(defparameter *step-height* 0.5)
(defparameter *channel-height* 1.0)
(defparameter *step-location* 1.0)
(defparameter *domain-length* 10.0)

(simulation "backward-step" :dim 2

  (domain (box (0 0) (*domain-length* *channel-height*))
    :mesh :n (200 20))

  (boundaries
    (inlet     (and (= x 0) (>= y *step-height*)))
    (outlet    (= x *domain-length*))
    (top       (= y *channel-height*))
    (bottom    (or (and (< x *step-location*) (= y *step-height*))
                   (and (>= x *step-location*) (= y 0))))
    (step-face (and (= x *step-location*) (< y *step-height*))))

  (physics :navier-stokes
    :formulation :velocity-pressure
    :Re 100

    (bc inlet :velocity
      :u (fn (x y)
           ;; Parabolic profile on inlet section
           (let* ((h (- *channel-height* *step-height*))
                  (y-local (- y *step-height*)))
             (* 1.5 (- 1 (expt (/ (- y-local (/ h 2)) (/ h 2)) 2)))))
      :v 0)
    (bc outlet    :outflow)
    (bc top       :no-slip)
    (bc bottom    :no-slip)
    (bc step-face :no-slip))

  (method :dcpse
    :neighbors 25
    :support-radius 5.0)

  (time :method :explicit-euler
    :dt 0.001
    :end 20.0)

  (output :vtk
    :directory "results/backward-step/"
    :every 1.0
    :fields (velocity pressure vorticity)))
