;;;; demo/st_cylinder.lisp - Space-Time Cylinder-in-Channel Flow
;;;;
;;;; StreamVorti SDL Example: Space-Time Navier-Stokes
;;;;
;;;; DFG benchmark (Schafer & Turek, 1996):
;;;;   Channel: 2.2 x 0.41
;;;;   Cylinder center: (0.2, 0.2), radius: 0.05
;;;;   Re = 20 (steady) or Re = 100 (unsteady / vortex shedding)
;;;;
;;;; Uses the ST-VMS formulation with SUPG/PSPG/LSIC stabilization.
;;;; For 2D spatial problems, the space-time domain is 3D with
;;;; tetrahedral (simplex) or hexahedral (tensor-product) elements.
;;;;
;;;; Run with:
;;;;   ./SpaceTimeNS -f demo/st_cylinder.lisp -lp lisp -pv
;;;;
;;;; Or generate mesh first:
;;;;   gmsh -2 demo/meshes/cylinder_in_channel.geo -format msh2 \
;;;;        -o demo/meshes/cylinder_in_channel.msh
;;;;   ./SpaceTimeNS -m demo/meshes/cylinder_in_channel.msh \
;;;;        -Re 100 -dt 0.1 -tf 8.0 -ov 2 -op 1 -et 0 -pv

(in-package :sdl)

;; Domain parameters (DFG benchmark)
(defparameter *channel-length* 2.2)
(defparameter *channel-height* 0.41)
(defparameter *cylinder-x* 0.2)
(defparameter *cylinder-y* 0.2)
(defparameter *cylinder-radius* 0.05)
(defparameter *U-max* 1.5)    ; Maximum inflow velocity

;; Parabolic inflow profile: u(y) = 4 * U_max * y * (H - y) / H^2
(defun parabolic-inflow (y)
  (* 4 *U-max* y (- *channel-height* y)
     (/ 1.0 (* *channel-height* *channel-height*))))

(simulation "cylinder-in-channel-st" :dim 2

  ;; Spatial domain: channel with cylindrical obstacle
  (domain (difference
            (box (0 0) (*channel-length* *channel-height*))
            (ball (*cylinder-x* *cylinder-y*) *cylinder-radius*))
    :mesh :generator :gmsh
    :file "demo/meshes/cylinder_in_channel.geo"
    :h 0.02)

  ;; Named boundaries
  (boundaries
    (inlet    (= x 0))
    (outlet   (= x *channel-length*))
    (top      (= y *channel-height*))
    (bottom   (= y 0))
    (cylinder (on-surface "cylinder")))

  ;; Space-time Navier-Stokes (velocity-pressure formulation)
  (physics :navier-stokes
    :formulation :velocity-pressure
    :method :space-time
    :Re 100

    ;; Boundary conditions
    (bc inlet    :inflow (lambda (x y t)
                           (values (parabolic-inflow y) 0)))
    (bc outlet   :outflow)
    (bc top      :no-slip)
    (bc bottom   :no-slip)
    (bc cylinder :no-slip))

  ;; Space-time discretization
  (spatial :fem
    :element-type :simplex         ; tet elements in space-time
    :velocity-order 2              ; P2 velocity
    :pressure-order 1)             ; P1 pressure

  ;; Space-time formulation options
  (space-time
    :formulation :st-vms           ; ST-VMS (Variational Multiscale)
    :temporal-scheme :dg           ; Discontinuous Galerkin in time
    :dt-slab 0.1                   ; Time slab thickness
    :time-steps-per-slab 1         ; Subdivisions within slab
    :stabilization (:supg :pspg :lsic))

  ;; Adaptive mesh refinement
  (adaptivity
    :strategy :space-time
    :indicator :residual           ; Element residual error indicator
    :refine-fraction 0.3
    :coarsen-fraction 0.05
    :max-level 5
    :max-elements 200000
    :frequency 2)                  ; Re-mesh every 2 slabs

  ;; Nonlinear solver (Newton)
  (nonlinear-solver
    :type :newton
    :max-iterations 30
    :rel-tolerance 1e-8
    :abs-tolerance 1e-12)

  ;; Linear solver
  (linear-solver
    :type :gmres
    :max-iterations 500
    :rel-tolerance 1e-10
    :preconditioner :ilu
    :gmres-restart 200)

  ;; VTK output for ParaView
  (output :vtk
    :directory "results/st_cylinder/"
    :every 1                       ; Save every time slab
    :fields (velocity pressure vorticity)))
