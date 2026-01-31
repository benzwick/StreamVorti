;;;; simulation.lisp - High-level simulation interface
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.

(in-package :sdl)

;;; ============================================================
;;; Simulation Accessors (for C++ interop)
;;; ============================================================

(defun simulation-name (sim)
  "Get simulation name."
  (simulation-data-name sim))

(defun simulation-dimension (sim)
  "Get simulation dimension."
  (simulation-data-dimension sim))

(defun simulation-mesh (sim)
  "Get simulation mesh specification."
  (simulation-data-mesh-spec sim))

(defun simulation-boundaries (sim)
  "Get simulation boundary conditions."
  (simulation-data-boundaries sim))

(defun simulation-physics (sim)
  "Get simulation physics parameters."
  (simulation-data-physics sim))

(defun simulation-discretization (sim)
  "Get simulation discretization parameters."
  (simulation-data-discretization sim))

(defun simulation-solver (sim)
  "Get simulation solver parameters."
  (simulation-data-solver sim))

(defun simulation-output (sim)
  "Get simulation output configuration."
  (simulation-data-output sim))

;;; ============================================================
;;; Simulation Validation
;;; ============================================================

(defun validate-simulation (sim)
  "Validate simulation configuration.
   Returns list of errors, empty if valid."
  (let ((errors nil))
    ;; Check required fields
    (unless (simulation-data-mesh-spec sim)
      (push "Missing mesh specification" errors))

    (unless (simulation-data-boundaries sim)
      (push "No boundary conditions defined" errors))

    ;; Check dimension consistency
    (let ((dim (simulation-data-dimension sim)))
      (unless (member dim '(2 3))
        (push (format nil "Invalid dimension: ~A (must be 2 or 3)" dim) errors)))

    (nreverse errors)))

;;; ============================================================
;;; Simulation Execution (Lisp-side, for future use)
;;; ============================================================

(defun realize-simulation (sim)
  "Convert simulation specification to runtime objects.
   Creates actual mesh and prepares for solving."
  (let ((mesh-spec (simulation-data-mesh-spec sim)))
    (when mesh-spec
      (setf (simulation-data-mesh sim)
            (streamvorti.mesh:realize-mesh-spec mesh-spec))))
  sim)

(defun print-simulation-summary (sim &optional (stream *standard-output*))
  "Print a summary of the simulation configuration."
  (format stream "~&Simulation: ~A~%" (simulation-data-name sim))
  (format stream "  Version: ~D~%" (simulation-data-version sim))
  (format stream "  Dimension: ~DD~%" (simulation-data-dimension sim))

  (let ((mesh-spec (simulation-data-mesh-spec sim)))
    (when mesh-spec
      (format stream "  Mesh: ~A~%"
              (streamvorti.mesh:mesh-spec-type mesh-spec))))

  (format stream "  Boundaries: ~D defined~%"
          (length (simulation-data-boundaries sim)))

  (let ((physics (simulation-data-physics sim)))
    (when physics
      (format stream "  Physics: ~A (~A)~%"
              (physics-data-type physics)
              (physics-data-formulation physics))
      (format stream "    Reynolds: ~A~%" (physics-data-reynolds physics))))

  (let ((disc (simulation-data-discretization sim)))
    (when disc
      (format stream "  Discretization: ~A~%"
              (discretization-data-method disc))
      (format stream "    Neighbors: ~D~%"
              (discretization-data-num-neighbors disc))))

  (let ((solver (simulation-data-solver sim)))
    (when solver
      (format stream "  Solver: ~A~%"
              (solver-data-timestepping solver))
      (format stream "    dt: ~A, end: ~A~%"
              (solver-data-dt solver)
              (solver-data-end-time solver))))

  sim)

;;; ============================================================
;;; File Loading Utilities
;;; ============================================================

(defun load-sdl-file (path)
  "Load an SDL file and return the simulation object."
  (setf *current-simulation* nil)
  (load path)
  (or *current-simulation*
      (error "No simulation defined in ~A" path)))

;;; ============================================================
;;; REPL Helpers
;;; ============================================================

(defun describe-current-simulation ()
  "Describe the current simulation in the REPL."
  (if *current-simulation*
      (print-simulation-summary *current-simulation*)
      (format t "~&No current simulation. Use (simulation ...) to define one.~%")))

(defun clear-simulation ()
  "Clear the current simulation."
  (setf *current-simulation* nil))
