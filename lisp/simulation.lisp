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
  spatial temporal output-config coupling)

(defvar *known-physics-types*
  '(:navier-stokes :heat-conduction :stokes :euler :laplace :diffusion
    :elasticity :poisson)
  "Known physics type keywords for parsing physics arguments.")

;;; ============================================================
;;; Compile-Time Walkers for Data-Oriented Simulation Macro
;;; ============================================================

(defun compile-bare-list (form)
  "Compile a bare list like (0 0) or (*len* *h*) into (list 0 0) etc.
   Atoms pass through unchanged."
  (if (atom form)
      form
      `(list ,@form)))

(defun compile-geom (form)
  "Compile a geometry expression from data syntax into constructor calls.
   (box (0 0) (1 1)) -> (box (list 0 0) (list 1 1))
   (ball (cx cy) r)  -> (ball (list cx cy) r)
   (difference g1 g2) -> (difference <compiled-g1> <compiled-g2>)"
  (cond
    ;; Atom: variable reference or literal, pass through
    ((atom form) form)
    ;; List with symbol car: dispatch on name
    ((symbolp (car form))
     (let ((name (string-upcase (symbol-name (car form)))))
       (cond
         ((string= name "BOX")
          `(box ,(compile-bare-list (second form))
                ,(compile-bare-list (third form))))
         ((string= name "BALL")
          `(ball ,(compile-bare-list (second form))
                 ,(third form)))
         ((string= name "DIFFERENCE")
          `(difference ,(compile-geom (second form))
                       ,(compile-geom (third form))))
         ((string= name "UNION")
          `(union-shape ,(compile-geom (second form))
                        ,(compile-geom (third form))))
         ((string= name "INTERSECTION")
          `(intersection-shape ,(compile-geom (second form))
                               ,(compile-geom (third form))))
         ;; Unknown symbol-headed list: treat as bare list
         (t `(list ,@form)))))
    ;; Bare list starting with a number
    ((numberp (car form))
     `(list ,@form))
    ;; Fallback
    (t form)))

(defun compile-domain-kwargs (kwargs)
  "Compile domain keyword arguments.
   :n (40 40) -> :n (list 40 40)
   Other values pass through."
  (let ((result nil))
    (loop with r = kwargs
          while r
          do (let ((key (car r))
                   (val (cadr r)))
               (push key result)
               (cond
                 ;; :n followed by a list → (list ...)
                 ((and (eq key :n) (listp val) (not (null val)))
                  (push `(list ,@val) result))
                 ;; Other values pass through
                 (t (push val result)))
               (setf r (cddr r))))
    (nreverse result)))

(defun compile-predicate (expr)
  "Compile a predicate expression into code that builds the predicate list
   at runtime.  Operators and coordinate symbols are quoted; value
   expressions are evaluated at runtime (allowing variable references).

   (= y 1)                     -> (list '= 'y 1)
   (and (= x 0) (>= y *h*))   -> (list 'and (list '= 'x 0) (list '>= 'y *h*))"
  (cond
    ((null expr) nil)
    ((numberp expr) expr)
    ((stringp expr) expr)
    ((symbolp expr)
     (let ((name (string-upcase (symbol-name expr))))
       (if (member name '("X" "Y" "Z" "=" ">" "<" ">=" "<=" "AND" "OR" "NOT")
                   :test #'string=)
           `',expr
           expr)))
    ((listp expr)
     `(list ,@(mapcar #'compile-predicate expr)))
    (t expr)))

(defun compile-comp (expr x y z)
  "Compile a single velocity component expression for use inside a lambda.
   Numbers pass through.  (fn (args) body) becomes an inline funcall.
   Other forms are funcalled with the coordinate variables."
  (cond
    ((null expr) 0)
    ((numberp expr) expr)
    ((and (listp expr)
          (symbolp (car expr))
          (string= (symbol-name (car expr)) "FN"))
     ;; (fn (params...) body...) → (funcall (lambda (params...) body...) x y)
     (let ((params (second expr)))
       `(funcall (lambda ,@(cdr expr))
                 ,x ,y ,@(when (> (length params) 2) (list z)))))
    (t `(funcall ,expr ,x ,y ,z))))

(defun compile-fn-expr (form)
  "Compile (fn (args) body...) into (lambda (args) body...).
   Non-fn forms pass through unchanged."
  (if (and (listp form)
           (symbolp (car form))
           (string= (symbol-name (car form)) "FN"))
      `(lambda ,@(cdr form))
      form))

(defun compile-bc-components (name type rest)
  "Compile a component-mode BC: :u expr :v expr [:w expr]
   into a single %make-bc with a combined lambda."
  (let ((u nil) (v nil) (w nil))
    (loop with r = rest
          while r
          do (let ((key (string-upcase (symbol-name (car r)))))
               (cond
                 ((string= key "U") (setf u (cadr r)) (setf r (cddr r)))
                 ((string= key "V") (setf v (cadr r)) (setf r (cddr r)))
                 ((string= key "W") (setf w (cadr r)) (setf r (cddr r)))
                 (t (setf r (cdr r))))))
    (let ((xv (gensym "X")) (yv (gensym "Y")) (zv (gensym "Z")))
      `(%make-bc :boundary ',name :type ,type
                 :function (lambda (,xv ,yv &optional ,zv)
                             ,@(unless w `((declare (ignore ,zv))))
                             (list ,(compile-comp u xv yv zv)
                                   ,(compile-comp v xv yv zv)
                                   ,@(when w
                                       (list (compile-comp w xv yv zv)))))))))

(defun compile-bc-convection (name type rest)
  "Compile a convection BC: :h val :T-inf val."
  (let ((h nil) (T-inf nil))
    (loop with r = rest
          while r
          do (cond
               ((eq (car r) :h)
                (setf h (cadr r)) (setf r (cddr r)))
               ((eq (car r) :T-inf)
                (setf T-inf (cadr r)) (setf r (cddr r)))
               (t (setf r (cdr r)))))
    `(%make-bc :boundary ',name :type ,type
               ,@(when h `(:h (coerce ,h 'double-float)))
               ,@(when T-inf `(:T-inf (coerce ,T-inf 'double-float))))))

(defun compile-bc (form)
  "Compile a (bc name type &rest args) form into a %make-bc call.

   Patterns:
   (bc wall :no-slip)                         -> type-only
   (bc outlet :pressure 0)                    -> scalar value
   (bc lid :velocity (1 0))                   -> bare list value
   (bc lid :velocity :u (fn ...) :v 0)        -> component functions
   (bc surface :convection :h 10 :T-inf 20)   -> keyword args
   (bc inlet :velocity :fn expr)              -> general function"
  (let* ((name (second form))
         (type (third form))
         (rest (cdddr form)))
    (cond
      ;; No args: type-only (no-slip, outflow, insulated, etc.)
      ((null rest)
       `(%make-bc :boundary ',name :type ,type))

      ;; Component mode: :u ... :v ... [:w ...]
      ((and (keywordp (car rest))
            (member (symbol-name (car rest)) '("U" "V" "W") :test #'string=))
       (compile-bc-components name type rest))

      ;; Convection mode: :h ... :T-inf ...
      ((and (keywordp (car rest))
            (string= (symbol-name (car rest)) "H"))
       (compile-bc-convection name type rest))

      ;; Function mode: :fn expr
      ((and (keywordp (car rest))
            (string= (symbol-name (car rest)) "FN"))
       `(%make-bc :boundary ',name :type ,type
                  :function ,(compile-fn-expr (second rest))))

      ;; Bare list value: (1 0) where car is a number
      ((and (listp (car rest))
            (numberp (caar rest)))
       `(%make-bc :boundary ',name :type ,type
                  :value (list ,@(car rest))))

      ;; Scalar value
      (t
       `(%make-bc :boundary ',name :type ,type
                  :value ,(car rest))))))

;;; --- Top-level form compilers ---

(defun compile-sim-domain (sim-var form)
  "Compile a (domain ...) form."
  (let ((args (cdr form)))
    (cond
      ;; File domain: (domain :file "path")
      ((and (keywordp (car args)) (eq (car args) :file))
       `(setf (sim-domain ,sim-var)
              (%make-domain :type :file :file ,(second args))))
      ;; Geometry domain: (domain geom-expr type &rest kwargs)
      (t
       (let* ((geom-form (car args))
              (rest (cdr args))
              (type (car rest))
              (kwargs (cdr rest)))
         `(setf (sim-domain ,sim-var)
                (domain ,(compile-geom geom-form) ,type
                        ,@(compile-domain-kwargs kwargs))))))))

(defun compile-sim-boundaries (sim-var form)
  "Compile a (boundaries (name pred) ...) form."
  (let* ((specs (cdr form))
         (compiled-specs
           (mapcar (lambda (spec)
                     `(list ',(car spec) ,(compile-predicate (cadr spec))))
                   specs)))
    `(setf (sim-boundaries ,sim-var)
           (make-boundaries (list ,@compiled-specs)))))

(defun compile-sim-subdomains (sim-var form)
  "Compile a (subdomains (name pred) ...) form."
  (let* ((specs (cdr form))
         (compiled-specs
           (mapcar (lambda (spec)
                     `(list ',(car spec) ,(compile-predicate (cadr spec))))
                   specs)))
    `(setf (sim-subdomains ,sim-var)
           (make-subdomains (list ,@compiled-specs)))))

(defun compile-sim-physics (sim-var form)
  "Compile a (physics [name] type {key val}* {(bc ...)}*) form."
  (let* ((args (cdr form))
         (rest args)
         (name nil) (type nil)
         (Re nil) (formulation nil) (conductivity nil) (subdomain nil)
         (bc-forms nil))
    ;; Parse name and type
    (cond
      ;; Two keywords, second is a known type -> first is name
      ((and (>= (length rest) 2)
            (keywordp (first rest))
            (keywordp (second rest))
            (member (second rest) *known-physics-types*))
       (setf name (first rest))
       (setf type (second rest))
       (setf rest (cddr rest)))
      ;; First is a keyword -> type
      ((and (>= (length rest) 1)
            (keywordp (first rest)))
       (setf type (first rest))
       (setf rest (cdr rest)))
      (t (error "Invalid physics form: ~S" form)))
    ;; Parse keyword args and (bc ...) forms
    (loop while rest
          do (cond
               ;; BC form: (bc ...)
               ((and (listp (car rest))
                     (symbolp (caar rest))
                     (string= (symbol-name (caar rest)) "BC"))
                (push (car rest) bc-forms)
                (setf rest (cdr rest)))
               ;; Keyword arg
               ((keywordp (car rest))
                (let ((key (car rest)) (val (cadr rest)))
                  (cond
                    ((eq key :Re) (setf Re val))
                    ((eq key :formulation) (setf formulation val))
                    ((eq key :conductivity) (setf conductivity val))
                    ((eq key :subdomain) (setf subdomain val)))
                  (setf rest (cddr rest))))
               (t (setf rest (cdr rest)))))
    (setf bc-forms (nreverse bc-forms))
    ;; Generate %make-physics call
    (let ((compiled-bcs (mapcar #'compile-bc bc-forms)))
      `(push (%make-physics
              ,@(when name `(:name ,name))
              :type ,type
              ,@(when Re `(:Re (coerce ,Re 'double-float)))
              ,@(when formulation `(:formulation ,formulation))
              ,@(when conductivity
                  `(:conductivity (coerce ,conductivity 'double-float)))
              ,@(when subdomain `(:subdomain ',subdomain))
              :bcs ,(if compiled-bcs
                        `(list ,@compiled-bcs)
                        nil))
             (sim-physics-list ,sim-var)))))

(defun compile-sim-spatial (sim-var form)
  "Compile a (spatial type &key ...) form."
  (let* ((args (cdr form))
         (type (car args))
         (kwargs (cdr args)))
    `(setf (sim-spatial ,sim-var) (make-spatial ,type ,@kwargs))))

(defun compile-sim-temporal (sim-var form)
  "Compile a (temporal type &key dt end tolerance) form.
   Type is the first positional argument, symmetric with spatial."
  (let* ((args (cdr form))
         (type (car args))
         (dt nil) (end nil) (tolerance nil))
    (loop with r = (cdr args)
          while r
          do (cond
               ((eq (car r) :dt)
                (setf dt (cadr r)) (setf r (cddr r)))
               ((eq (car r) :end)
                (setf end (cadr r)) (setf r (cddr r)))
               ((eq (car r) :tolerance)
                (setf tolerance (cadr r)) (setf r (cddr r)))
               (t (setf r (cdr r)))))
    `(setf (sim-temporal ,sim-var)
           (%make-temporal
            :method ,type
            ,@(when dt `(:dt (coerce ,dt 'double-float)))
            ,@(when end `(:end (coerce ,end 'double-float)))
            ,@(when tolerance
                `(:tolerance (coerce ,tolerance 'double-float)))))))

(defun compile-sim-output (sim-var form)
  "Compile a (output format &key directory every fields probes) form."
  (let ((fmt (second form))
        (directory nil) (every nil) (fields nil) (probes nil))
    (loop with r = (cddr form)
          while r
          do (cond
               ((eq (car r) :directory)
                (setf directory (cadr r)) (setf r (cddr r)))
               ((eq (car r) :every)
                (setf every (cadr r)) (setf r (cddr r)))
               ((eq (car r) :fields)
                (setf fields (cadr r)) (setf r (cddr r)))
               ((eq (car r) :probes)
                (setf probes (cadr r)) (setf r (cddr r)))
               (t (setf r (cdr r)))))
    `(setf (sim-output-config ,sim-var)
           (%make-output
            :format ,fmt
            ,@(when directory `(:directory ,directory))
            ,@(when every `(:interval (coerce ,every 'double-float)))
            ,@(when fields
                `(:fields ',(if (listp fields) fields (list fields))))
            ,@(when probes `(:probes ',probes))))))

(defun compile-sim-coupling (sim-var form)
  "Compile a (coupling name &key physics type iterations interface) form."
  (let ((name (second form))
        (physics nil) (type nil) (iterations nil) (interface nil))
    (loop with r = (cddr form)
          while r
          do (cond
               ((eq (car r) :physics)
                (setf physics (cadr r)) (setf r (cddr r)))
               ((eq (car r) :type)
                (setf type (cadr r)) (setf r (cddr r)))
               ((eq (car r) :iterations)
                (setf iterations (cadr r)) (setf r (cddr r)))
               ((eq (car r) :interface)
                (setf interface (cadr r)) (setf r (cddr r)))
               (t (setf r (cdr r)))))
    `(setf (sim-coupling ,sim-var)
           (%make-coupling
            :name ,name
            ,@(when physics
                `(:physics ',(if (listp physics) physics (list physics))))
            ,@(when type `(:type ,type))
            ,@(when iterations `(:iterations ,iterations))
            ,@(when interface `(:interface ',interface))))))

(defun compile-sim-form (sim-var form)
  "Compile a single simulation body form by dispatching on (car form)."
  (unless (and (listp form) (symbolp (car form)))
    (error "Invalid simulation body form: ~S" form))
  (let ((tag (string-upcase (symbol-name (car form)))))
    (cond
      ((string= tag "DOMAIN")      (compile-sim-domain sim-var form))
      ((string= tag "BOUNDARIES")  (compile-sim-boundaries sim-var form))
      ((string= tag "SUBDOMAINS")  (compile-sim-subdomains sim-var form))
      ((string= tag "PHYSICS")     (compile-sim-physics sim-var form))
      ((string= tag "SPATIAL")     (compile-sim-spatial sim-var form))
      ((string= tag "TEMPORAL")    (compile-sim-temporal sim-var form))
      ((string= tag "OUTPUT")      (compile-sim-output sim-var form))
      ((string= tag "COUPLING")    (compile-sim-coupling sim-var form))
      (t (error "Unknown simulation form: ~A" (car form))))))

;;; ============================================================
;;; Simulation Macro
;;; ============================================================

(defmacro simulation (name &rest args)
  "Define a complete simulation using data-oriented syntax.

   The body forms are walked at compile time and transformed into
   constructor calls.  Symbols are matched by SYMBOL-NAME, so the
   macro works regardless of the caller's package.

   Example:
   (simulation \"lid-driven-cavity\" :dim 2
     (domain (box (0 0) (1 1)) :mesh :n (40 40))
     (boundaries
       (lid    (= y 1))
       (bottom (= y 0))
       (left   (= x 0))
       (right  (= x 1)))
     (physics :navier-stokes :Re 100
       (bc lid    :velocity (1 0))
       (bc bottom :no-slip)
       (bc left   :no-slip)
       (bc right  :no-slip))
     (temporal :explicit-euler :dt 0.001 :end 10.0)
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
         ,@(mapcar (lambda (form) (compile-sim-form sim-var form))
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

(defun simulation-spatial (sim)
  "Get simulation spatial discretization."
  (sim-spatial sim))

(defun simulation-temporal (sim)
  "Get simulation temporal integration."
  (sim-temporal sim))

(defun simulation-coupling (sim)
  "Get simulation coupling definition."
  (sim-coupling sim))
