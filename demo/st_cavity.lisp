;;;; demo/st_cavity.lisp - Space-Time Lid-Driven Cavity Flow
;;;;
;;;; StreamVorti SDL Example: Space-Time Navier-Stokes
;;;;
;;;; Classic lid-driven cavity using space-time FEM.
;;;; Reference: Ghia, Ghia & Shin (1982)
;;;;
;;;; Run with:
;;;;   ./SpaceTimeNS -f demo/st_cavity.lisp -lp lisp -pv
;;;;
;;;; Or CLI mode:
;;;;   ./SpaceTimeNS -m demo/meshes/cavity_quad.mesh -Re 100 \
;;;;        -dt 0.5 -tf 30.0 -ov 2 -op 1 -et 1 -rs 2 -pv

(in-package :sdl)

(simulation "lid-driven-cavity-st" :dim 2

  ;; Unit square domain with structured quad mesh
  (domain (box (0 0) (1 1))
    :mesh :n (20 20)
    :element-type :quad)           ; Quads -> hexahedra in space-time

  ;; Named boundaries
  (boundaries
    (bottom (= y 0))
    (right  (= x 1))
    (lid    (= y 1))
    (left   (= x 0)))

  ;; Space-time Navier-Stokes
  (physics :navier-stokes
    :formulation :velocity-pressure
    :method :space-time
    :Re 100

    (bc bottom :no-slip)
    (bc right  :no-slip)
    (bc lid    :velocity (1 0))
    (bc left   :no-slip))

  ;; Finite element discretization
  (spatial :fem
    :element-type :tensor-product  ; hex elements in space-time
    :velocity-order 2
    :pressure-order 1)

  ;; Space-time options
  (space-time
    :formulation :st-vms
    :temporal-scheme :dg
    :dt-slab 0.5
    :time-steps-per-slab 1
    :stabilization (:supg :pspg :lsic))

  ;; Adaptive refinement
  (adaptivity
    :strategy :spatial             ; Refine only in space
    :indicator :vorticity          ; Refine near vortex cores
    :refine-fraction 0.2
    :max-level 4
    :max-elements 100000
    :frequency 5)

  ;; Newton solver
  (nonlinear-solver
    :type :newton
    :max-iterations 20
    :rel-tolerance 1e-8)

  ;; Direct solver (good for small problems)
  (linear-solver
    :type :umfpack)

  ;; Line probes for Ghia validation
  (probes
    (line vertical-center   (x 0.5))
    (line horizontal-center (y 0.5)))

  ;; Output
  (output :vtk
    :directory "results/st_cavity/"
    :every 1
    :fields (velocity pressure vorticity)))
