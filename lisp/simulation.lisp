;;;; simulation.lisp - SDL v2 simulation assembly
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
;;; Simulation Structure
;;; ============================================================

(defstruct (sim-data (:constructor %make-sim)
                     (:conc-name sim-))
  name dim domain boundaries subdomains physics-list
  time-solver output-config coupling)

;;; ============================================================
;;; Simulation Helper Functions
;;; ============================================================

(defun boundaries (specs)
  "Create a boundary set for use in a simulation definition.
   Returns a tagged wrapper so simulation can identify it.

   Example:
   (boundaries '((lid (= y 1))
                 (walls (or (= x 0) (= x 1)))))"
  (%make-boundary-set :defs (make-boundaries specs)))

(defvar *known-physics-types*
  '(:navier-stokes :heat-conduction :stokes :euler :laplace :diffusion
    :elasticity :poisson)
  "Known physics type keywords for parsing physics helper arguments.")

(defun physics (&rest args)
  "Create a physics definition for use in a simulation definition.

   Single physics:
   (physics :navier-stokes :Re 100 '((lid :velocity (1 0)) ...))

   Named physics (for multiphysics):
   (physics :flow :navier-stokes :Re 100 :subdomain 'fluid '((lid ...) ...))
   (physics :thermal :heat-conduction '((inlet :temperature 20) ...))"
  (let* (;; Last arg is always the BCs list
         (bc-specs (car (last args)))
         (rest-args (butlast args))
         ;; Parse name and type
         (name nil) (type nil) (keywords nil))
    ;; Determine if first arg is a name or a type
    (cond
      ;; Two keywords at start, second is a known type -> first is name
      ((and (>= (length rest-args) 2)
            (keywordp (first rest-args))
            (keywordp (second rest-args))
            (member (second rest-args) *known-physics-types*))
       (setf name (first rest-args))
       (setf type (second rest-args))
       (setf keywords (cddr rest-args)))
      ;; First is a keyword (known or not) -> it's the type
      ((and (>= (length rest-args) 1)
            (keywordp (first rest-args)))
       (setf type (first rest-args))
       (setf keywords (cdr rest-args)))
      (t (error "Invalid physics arguments: ~S" args)))
    ;; Parse keyword pairs from remaining args
    (let ((Re nil) (conductivity nil) (subdomain nil))
      (loop with r = keywords
            while r
            do (cond
                 ((eq (car r) :Re)
                  (setf Re (cadr r)) (setf r (cddr r)))
                 ((eq (car r) :conductivity)
                  (setf conductivity (cadr r)) (setf r (cddr r)))
                 ((eq (car r) :subdomain)
                  (setf subdomain (cadr r)) (setf r (cddr r)))
                 (t (setf r (cdr r)))))
      ;; Create physics with parsed BCs
      (%make-physics
       :name name
       :type type
       :Re (when Re (coerce Re 'double-float))
       :conductivity (when conductivity (coerce conductivity 'double-float))
       :subdomain subdomain
       :bcs (when bc-specs (mapcar #'parse-bc-spec bc-specs))))))

(defun time (&key dt end (method nil))
  "Create a time solver for use in a simulation definition.

   Example:
   (time :dt 0.001 :end 10.0)"
  (%make-solver
   :method method
   :dt (when dt (coerce dt 'double-float))
   :end (when end (coerce end 'double-float))))

(defun output (format &rest args)
  "Create output configuration for use in a simulation definition.

   Example:
   (output :vtk :every 0.1)"
  (let ((every nil) (directory nil) (fields nil) (probes nil))
    (loop with rest = args
          while rest
          do (cond
               ((eq (car rest) :every)
                (setf every (cadr rest)) (setf rest (cddr rest)))
               ((eq (car rest) :directory)
                (setf directory (cadr rest)) (setf rest (cddr rest)))
               ((eq (car rest) :fields)
                (setf fields (cadr rest)) (setf rest (cddr rest)))
               ((eq (car rest) :probes)
                (setf probes (cadr rest)) (setf rest (cddr rest)))
               (t (setf rest (cdr rest)))))
    (%make-output
     :format format
     :directory directory
     :interval (when every (coerce every 'double-float))
     :fields fields
     :probes probes)))

(defun subdomains (specs)
  "Create a subdomain set for use in a simulation definition.

   Example:
   (subdomains '((fluid (< x 0.5))
                 (solid (>= x 0.5))))"
  (%make-subdomain-set :defs (make-subdomains specs)))

(defun coupling (name &key physics type iterations interface)
  "Create a coupling definition for use in a simulation definition.

   Example:
   (coupling :flow-thermal :physics '(flow thermal) :type :sequential)"
  (%make-coupling
   :name name
   :physics physics
   :type type
   :iterations iterations
   :interface interface))

;;; ============================================================
;;; Simulation Component Dispatch
;;; ============================================================

(defun %sim-add (sim obj)
  "Add a component to a simulation by dispatching on type."
  (typecase obj
    (domain-data
     (setf (sim-domain sim) obj))
    (boundary-set
     (setf (sim-boundaries sim) (boundary-set-defs obj)))
    (subdomain-set
     (setf (sim-subdomains sim) (subdomain-set-defs obj)))
    (physics-data
     (push obj (sim-physics-list sim)))
    (solver-data
     (setf (sim-time-solver sim) obj))
    (output-data
     (setf (sim-output-config sim) obj))
    (coupling-data
     (setf (sim-coupling sim) obj))
    (t (warn "Unknown simulation component type: ~A" (type-of obj)))))

;;; ============================================================
;;; Simulation Macro
;;; ============================================================

(defmacro simulation (name &rest args)
  "Define a complete simulation.

   Example:
   (simulation \"lid-driven-cavity\" :dim 2
     (domain (box '(0 0) '(1 1)) :mesh :n '(40 40))
     (boundaries '((lid (= y 1)) (walls (or (= x 0) (= x 1)))))
     (physics :navier-stokes :Re 100 '((lid :velocity (1 0))))
     (time :dt 0.001 :end 10.0)
     (output :vtk :every 0.1))"
  (let ((dim 2)
        (body-forms nil))
    ;; Parse :dim from args, collect body forms
    (loop with rest = args
          while rest
          do (cond
               ((eq (car rest) :dim)
                (setf dim (cadr rest))
                (setf rest (cddr rest)))
               (t
                (push (car rest) body-forms)
                (setf rest (cdr rest)))))
    (setf body-forms (nreverse body-forms))
    (let ((sim-var (gensym "SIM")))
      `(let ((,sim-var (%make-sim :name ,name :dim ,dim)))
         ,@(mapcar (lambda (form)
                     `(%sim-add ,sim-var ,form))
                   body-forms)
         ;; Reverse physics list to preserve definition order
         (setf (sim-physics-list ,sim-var)
               (nreverse (sim-physics-list ,sim-var)))
         ,sim-var))))

;;; ============================================================
;;; Simulation Accessors
;;; ============================================================

(defun simulation-name (sim)
  "Get simulation name."
  (sim-name sim))

(defun simulation-dim (sim)
  "Get simulation dimension."
  (sim-dim sim))

(defun simulation-domain (sim)
  "Get simulation domain."
  (sim-domain sim))

(defun simulation-boundaries (sim)
  "Get simulation boundary definitions."
  (sim-boundaries sim))

(defun simulation-physics (sim)
  "Get the primary (first) physics definition."
  (first (sim-physics-list sim)))

(defun simulation-subdomains (sim)
  "Get simulation subdomain definitions."
  (sim-subdomains sim))

(defun simulation-physics-list (sim)
  "Get all physics definitions."
  (sim-physics-list sim))

(defun simulation-coupling (sim)
  "Get simulation coupling definition."
  (sim-coupling sim))
