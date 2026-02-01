;;;; test-sdl.lisp - Test suite for new SDL design (ADR-0001)
;;;;
;;;; Tests for the Simulation Definition Language redesign.
;;;; Run with:
;;;;   sbcl --load lisp/sdl.lisp --load lisp/test-sdl.lisp
;;;;   ecl -load lisp/sdl.lisp -load lisp/test-sdl.lisp

(defpackage :streamvorti.sdl.tests
  (:use :cl)
  (:export #:run-all-tests
           #:*pass-count*
           #:*fail-count*))

(in-package :streamvorti.sdl.tests)

;;; ============================================================
;;; Test Framework (minimal, no dependencies)
;;; ============================================================

(defvar *test-count* 0)
(defvar *pass-count* 0)
(defvar *fail-count* 0)
(defvar *current-test* nil)

(defmacro define-test (name &body body)
  "Define a test case."
  `(defun ,name ()
     (let ((*current-test* ',name))
       (format t "~%Testing ~A...~%" ',name)
       (handler-case
           (progn ,@body t)
         (error (e)
           (incf *fail-count*)
           (format t "  ERROR: ~A~%" e)
           nil)))))

(defun check (condition &optional (message ""))
  "Assert that CONDITION is true."
  (incf *test-count*)
  (if condition
      (progn
        (incf *pass-count*)
        t)
      (progn
        (incf *fail-count*)
        (format t "  FAIL: ~A~%" message)
        nil)))

(defun check-equal (actual expected &optional (message ""))
  "Assert that ACTUAL equals EXPECTED."
  (check (equal actual expected)
         (format nil "~A~%    Expected: ~S~%    Actual:   ~S"
                 message expected actual)))

(defun check-near (actual expected &optional (tolerance 1e-6) (message ""))
  "Assert that ACTUAL is within TOLERANCE of EXPECTED."
  (check (< (abs (- actual expected)) tolerance)
         (format nil "~A~%    Expected: ~F (±~F)~%    Actual:   ~F"
                 message expected tolerance actual)))

(defun check-type-is (obj type &optional (message ""))
  "Assert that OBJ is of TYPE."
  (check (typep obj type)
         (format nil "~A~%    Expected type: ~S~%    Actual type:   ~S"
                 message type (type-of obj))))

(defun check-not-null (obj &optional (message ""))
  "Assert that OBJ is not NIL."
  (check (not (null obj))
         (format nil "~A~%    Expected non-NIL, got NIL" message)))

(defun check-signals (condition-type thunk &optional (message ""))
  "Assert that THUNK signals a condition of CONDITION-TYPE."
  (incf *test-count*)
  (handler-case
      (progn
        (funcall thunk)
        (incf *fail-count*)
        (format t "  FAIL: ~A~%    Expected ~S to be signaled~%"
                message condition-type)
        nil)
    (condition (c)
      (if (typep c condition-type)
          (progn (incf *pass-count*) t)
          (progn
            (incf *fail-count*)
            (format t "  FAIL: ~A~%    Expected ~S, got ~S~%"
                    message condition-type (type-of c))
            nil)))))

;;; ============================================================
;;; Tests for Geometry Primitives
;;; ============================================================

(define-test test-box-creation
  "Test creating box geometry in various dimensions."
  ;; 1D box (interval)
  (let ((b1 (sdl:box '(0) '(1))))
    (check-not-null b1 "1D box created")
    (check-equal (sdl:geometry-dimension b1) 1 "1D box dimension"))

  ;; 2D box (rectangle)
  (let ((b2 (sdl:box '(0 0) '(1 1))))
    (check-not-null b2 "2D box created")
    (check-equal (sdl:geometry-dimension b2) 2 "2D box dimension")
    (check-equal (sdl:geometry-min b2) '(0 0) "2D box min")
    (check-equal (sdl:geometry-max b2) '(1 1) "2D box max"))

  ;; 3D box (cuboid)
  (let ((b3 (sdl:box '(0 0 0) '(1 2 3))))
    (check-not-null b3 "3D box created")
    (check-equal (sdl:geometry-dimension b3) 3 "3D box dimension")
    (check-equal (sdl:geometry-max b3) '(1 2 3) "3D box max")))

(define-test test-ball-creation
  "Test creating ball geometry in various dimensions."
  ;; 2D ball (circle)
  (let ((c (sdl:ball '(0.5 0.5) 0.25)))
    (check-not-null c "2D ball created")
    (check-equal (sdl:geometry-dimension c) 2 "2D ball dimension")
    (check-equal (sdl:ball-center c) '(0.5 0.5) "2D ball center")
    (check-near (sdl:ball-radius c) 0.25 1e-10 "2D ball radius"))

  ;; 3D ball (sphere)
  (let ((s (sdl:ball '(0 0 0) 1.0)))
    (check-not-null s "3D ball created")
    (check-equal (sdl:geometry-dimension s) 3 "3D ball dimension")))

(define-test test-geometry-contains
  "Test point containment for geometries."
  (let ((b (sdl:box '(0 0) '(1 1))))
    (check (sdl:geometry-contains-p b '(0.5 0.5)) "Center point in box")
    (check (sdl:geometry-contains-p b '(0 0)) "Corner point in box")
    (check (not (sdl:geometry-contains-p b '(1.5 0.5))) "Point outside box"))

  (let ((c (sdl:ball '(0 0) 1.0)))
    (check (sdl:geometry-contains-p c '(0 0)) "Center in ball")
    (check (sdl:geometry-contains-p c '(0.5 0.5)) "Point inside ball")
    (check (not (sdl:geometry-contains-p c '(1 1))) "Point outside ball")))

;;; ============================================================
;;; Tests for Domain Creation
;;; ============================================================

(define-test test-domain-mesh-structured
  "Test creating structured mesh domains."
  ;; 2D structured mesh on box
  (let ((d (sdl:domain (sdl:box '(0 0) '(1 1)) :mesh :n '(10 10))))
    (check-not-null d "Domain created")
    (check-equal (sdl:domain-type d) :mesh "Domain type is mesh")
    (check-equal (sdl:domain-dimension d) 2 "Domain dimension")
    (check-equal (sdl:domain-n d) '(10 10) "Domain divisions"))

  ;; 3D structured mesh
  (let ((d (sdl:domain (sdl:box '(0 0 0) '(1 1 1)) :mesh :n '(5 5 5))))
    (check-equal (sdl:domain-dimension d) 3 "3D domain dimension")
    (check-equal (sdl:domain-n d) '(5 5 5) "3D domain divisions")))

(define-test test-domain-points
  "Test creating point cloud domains."
  ;; Points with divisions
  (let ((d (sdl:domain (sdl:box '(0 0) '(1 1)) :points :n '(20 20))))
    (check-equal (sdl:domain-type d) :points "Domain type is points")
    (check-equal (sdl:domain-n d) '(20 20) "Point divisions"))

  ;; Points with spacing
  (let ((d (sdl:domain (sdl:box '(0 0) '(1 1)) :points :h 0.05)))
    (check-equal (sdl:domain-type d) :points "Domain type is points")
    (check-near (sdl:domain-h d) 0.05 1e-10 "Point spacing")))

(define-test test-domain-from-file
  "Test creating domain from file."
  (let ((d (sdl:domain :file "test.mesh")))
    (check-equal (sdl:domain-type d) :file "Domain type is file")
    (check-equal (sdl:domain-file d) "test.mesh" "Domain file path")))

(define-test test-domain-particles
  "Test creating particle domains."
  (let ((d (sdl:domain (sdl:box '(0 0) '(1 1)) :particles :n 1000)))
    (check-equal (sdl:domain-type d) :particles "Domain type is particles")
    (check-equal (sdl:domain-particle-count d) 1000 "Particle count")))

;;; ============================================================
;;; Tests for Boundary Definitions
;;; ============================================================

(define-test test-boundary-predicate-equality
  "Test boundary predicates for coordinate equality."
  ;; Simple equality predicate
  (let ((pred (sdl:make-predicate '(= y 1))))
    (check-not-null pred "Predicate created")
    (check (sdl:predicate-matches-p pred '(0.5 1.0 0)) "Point on y=1 matches")
    (check (not (sdl:predicate-matches-p pred '(0.5 0.5 0))) "Point not on y=1"))

  ;; x equality
  (let ((pred (sdl:make-predicate '(= x 0))))
    (check (sdl:predicate-matches-p pred '(0 0.5 0)) "Point on x=0 matches")
    (check (not (sdl:predicate-matches-p pred '(0.5 0.5 0))) "Point not on x=0")))

(define-test test-boundary-predicate-boolean
  "Test boolean combinations of predicates."
  ;; OR predicate
  (let ((pred (sdl:make-predicate '(or (= x 0) (= x 1)))))
    (check (sdl:predicate-matches-p pred '(0 0.5 0)) "x=0 matches")
    (check (sdl:predicate-matches-p pred '(1 0.5 0)) "x=1 matches")
    (check (not (sdl:predicate-matches-p pred '(0.5 0.5 0))) "x=0.5 doesn't match"))

  ;; AND predicate
  (let ((pred (sdl:make-predicate '(and (> x 0) (< x 1)))))
    (check (sdl:predicate-matches-p pred '(0.5 0.5 0)) "0<x<1 matches")
    (check (not (sdl:predicate-matches-p pred '(0 0.5 0))) "x=0 doesn't match"))

  ;; NOT predicate
  (let ((pred (sdl:make-predicate '(not (= y 1)))))
    (check (sdl:predicate-matches-p pred '(0.5 0.5 0)) "y!=1 matches")
    (check (not (sdl:predicate-matches-p pred '(0.5 1 0))) "y=1 doesn't match")))

(define-test test-boundary-attribute
  "Test boundary selection by mesh attribute."
  (let ((sel (sdl:make-boundary-selector '(attribute 3))))
    (check-equal (sdl:selector-type sel) :attribute "Selector type")
    (check-equal (sdl:selector-attribute sel) 3 "Attribute value")))

(define-test test-boundaries-block
  "Test the boundaries block syntax."
  (let ((boundaries (sdl:make-boundaries
                      '((lid    (= y 1))
                        (bottom (= y 0))
                        (walls  (or (= x 0) (= x 1)))))))
    (check-equal (length boundaries) 3 "Three boundaries defined")
    (check-equal (sdl:boundary-name (first boundaries)) 'lid "First boundary name")
    (check-equal (sdl:boundary-name (second boundaries)) 'bottom "Second boundary name")
    (check-equal (sdl:boundary-name (third boundaries)) 'walls "Third boundary name")))

;;; ============================================================
;;; Tests for Subdomains
;;; ============================================================

(define-test test-subdomains-block
  "Test the subdomains block syntax."
  (let ((subdomains (sdl:make-subdomains
                      '((fluid (< x 0.5))
                        (solid (>= x 0.5))))))
    (check-equal (length subdomains) 2 "Two subdomains defined")
    (check-equal (sdl:subdomain-name (first subdomains)) 'fluid "First subdomain name")))

;;; ============================================================
;;; Tests for Physics
;;; ============================================================

(define-test test-physics-navier-stokes
  "Test Navier-Stokes physics definition."
  (let ((phys (sdl:make-physics
                :type :navier-stokes
                :Re 100
                :bcs '((lid :velocity (1 0))
                       (walls :no-slip)))))
    (check-equal (sdl:physics-type phys) :navier-stokes "Physics type")
    (check-near (sdl:physics-Re phys) 100.0 1e-10 "Reynolds number")
    (check-equal (length (sdl:physics-bcs phys)) 2 "Two BCs defined")))

(define-test test-physics-heat-conduction
  "Test heat conduction physics definition."
  (let ((phys (sdl:make-physics
                :type :heat-conduction
                :conductivity 1.0
                :bcs '((lid :temperature 100)
                       (walls :insulated)))))
    (check-equal (sdl:physics-type phys) :heat-conduction "Physics type")
    (check-near (sdl:physics-conductivity phys) 1.0 1e-10 "Conductivity")))

(define-test test-boundary-conditions
  "Test various boundary condition types."
  ;; Velocity BC
  (let ((bc (sdl:make-bc 'inlet :velocity '(1 0 0))))
    (check-equal (sdl:bc-type bc) :velocity "BC type")
    (check-equal (sdl:bc-boundary bc) 'inlet "BC boundary name")
    (check-equal (sdl:bc-value bc) '(1 0 0) "BC value"))

  ;; No-slip BC
  (let ((bc (sdl:make-bc 'wall :no-slip)))
    (check-equal (sdl:bc-type bc) :no-slip "No-slip BC type"))

  ;; Pressure BC
  (let ((bc (sdl:make-bc 'outlet :pressure 0)))
    (check-equal (sdl:bc-type bc) :pressure "Pressure BC type")
    (check-equal (sdl:bc-value bc) 0 "Pressure value"))

  ;; Temperature BC
  (let ((bc (sdl:make-bc 'hot :temperature 100)))
    (check-equal (sdl:bc-type bc) :temperature "Temperature BC type")
    (check-equal (sdl:bc-value bc) 100 "Temperature value"))

  ;; Convection BC
  (let ((bc (sdl:make-bc 'surface :convection :h 10 :T-inf 20)))
    (check-equal (sdl:bc-type bc) :convection "Convection BC type")
    (check-near (sdl:bc-h bc) 10.0 1e-10 "Heat transfer coeff")
    (check-near (sdl:bc-T-inf bc) 20.0 1e-10 "Ambient temperature")))

(define-test test-bc-with-function
  "Test boundary conditions with user functions."
  (let ((bc (sdl:make-bc 'inlet :velocity
                         :fn (lambda (x y z) (list (* 4 y (- 1 y)) 0 0)))))
    (check-equal (sdl:bc-type bc) :velocity "Velocity BC type")
    (check-not-null (sdl:bc-function bc) "BC has function")
    ;; Evaluate at y=0.5 (parabolic max)
    (let ((v (funcall (sdl:bc-function bc) 0 0.5 0)))
      (check-near (first v) 1.0 1e-10 "Parabolic profile at y=0.5"))))

;;; ============================================================
;;; Tests for Methods (Discretization)
;;; ============================================================

(define-test test-method-dcpse
  "Test DCPSE method definition."
  (let ((m (sdl:make-method :dcpse :neighbors 25 :support-radius 3.5)))
    (check-equal (sdl:method-type m) :dcpse "Method type")
    (check-equal (sdl:method-neighbors m) 25 "Neighbor count")
    (check-near (sdl:method-support-radius m) 3.5 1e-10 "Support radius")))

(define-test test-method-fem
  "Test FEM method definition."
  (let ((m (sdl:make-method :fem :order 2)))
    (check-equal (sdl:method-type m) :fem "Method type")
    (check-equal (sdl:method-order m) 2 "Element order")))

(define-test test-method-sph
  "Test SPH method definition."
  (let ((m (sdl:make-method :sph :kernel :wendland :h 0.02)))
    (check-equal (sdl:method-type m) :sph "Method type")
    (check-equal (sdl:method-kernel m) :wendland "Kernel type")
    (check-near (sdl:method-h m) 0.02 1e-10 "Smoothing length")))

;;; ============================================================
;;; Tests for Solvers
;;; ============================================================

(define-test test-time-solver
  "Test time integration solver."
  (let ((s (sdl:make-time-solver :method :bdf2 :dt 0.001 :end 10.0)))
    (check-equal (sdl:solver-method s) :bdf2 "Time method")
    (check-near (sdl:solver-dt s) 0.001 1e-10 "Time step")
    (check-near (sdl:solver-end s) 10.0 1e-10 "End time")))

(define-test test-steady-solver
  "Test steady-state solver."
  (let ((s (sdl:make-steady-solver :method :newton :tolerance 1e-8)))
    (check-equal (sdl:solver-method s) :newton "Steady method")
    (check-near (sdl:solver-tolerance s) 1e-8 1e-15 "Tolerance")))

(define-test test-linear-solver
  "Test linear solver configuration."
  (let ((s (sdl:make-linear-solver :method :gmres :preconditioner :ilu)))
    (check-equal (sdl:linear-solver-method s) :gmres "Linear solver method")
    (check-equal (sdl:linear-solver-preconditioner s) :ilu "Preconditioner")))

;;; ============================================================
;;; Tests for Materials
;;; ============================================================

(define-test test-material-definition
  "Test material property definition."
  (let ((m (sdl:make-material 'steel
                              :density 7850
                              :youngs-modulus 200e9
                              :poissons-ratio 0.3)))
    (check-equal (sdl:material-name m) 'steel "Material name")
    (check-near (sdl:material-density m) 7850.0 1e-10 "Density")
    (check-near (sdl:material-youngs-modulus m) 200e9 1e-10 "Young's modulus")
    (check-near (sdl:material-poissons-ratio m) 0.3 1e-10 "Poisson's ratio")))

(define-test test-fluid-material
  "Test fluid material properties."
  (let ((m (sdl:make-material 'water
                              :density 1000
                              :viscosity 0.001
                              :conductivity 0.6)))
    (check-near (sdl:material-viscosity m) 0.001 1e-10 "Viscosity")
    (check-near (sdl:material-conductivity m) 0.6 1e-10 "Conductivity")))

;;; ============================================================
;;; Tests for Output Configuration
;;; ============================================================

(define-test test-output-basic
  "Test basic output configuration."
  (let ((o (sdl:make-output :format :vtk
                            :directory "results/"
                            :every 0.1
                            :fields '(velocity pressure))))
    (check-equal (sdl:output-format o) :vtk "Output format")
    (check-equal (sdl:output-directory o) "results/" "Output directory")
    (check-near (sdl:output-interval o) 0.1 1e-10 "Output interval")
    (check-equal (sdl:output-fields o) '(velocity pressure) "Output fields")))

(define-test test-output-probes
  "Test probe output configuration."
  (let ((o (sdl:make-output :format :vtk
                            :probes '((:point (0.5 0.5) :fields (velocity pressure))
                                      (:line (0 0.5) (1 0.5) :n 100 :fields (velocity))))))
    (check-equal (length (sdl:output-probes o)) 2 "Two probes defined")))

;;; ============================================================
;;; Tests for Coupling (Multiphysics)
;;; ============================================================

(define-test test-coupling-sequential
  "Test sequential multiphysics coupling."
  (let ((c (sdl:make-coupling :name :thermal-flow
                              :physics '(thermal flow)
                              :type :sequential
                              :iterations 10)))
    (check-equal (sdl:coupling-name c) :thermal-flow "Coupling name")
    (check-equal (sdl:coupling-type c) :sequential "Coupling type")
    (check-equal (sdl:coupling-physics c) '(thermal flow) "Coupled physics")
    (check-equal (sdl:coupling-iterations c) 10 "Coupling iterations")))

(define-test test-coupling-monolithic
  "Test monolithic multiphysics coupling."
  (let ((c (sdl:make-coupling :name :fsi
                              :physics '(fluid structure)
                              :type :monolithic
                              :interface 'moving-wall)))
    (check-equal (sdl:coupling-type c) :monolithic "Coupling type")
    (check-equal (sdl:coupling-interface c) 'moving-wall "Interface")))

;;; ============================================================
;;; Tests for Complete Simulation
;;; ============================================================

(define-test test-simulation-lid-driven-cavity
  "Test complete lid-driven cavity simulation definition."
  (let ((sim (sdl:simulation "lid-driven-cavity" :dim 2
               (sdl:domain (sdl:box '(0 0) '(1 1)) :mesh :n '(40 40))

               (sdl:boundaries
                 '((lid    (= y 1))
                   (bottom (= y 0))
                   (left   (= x 0))
                   (right  (= x 1))))

               (sdl:physics :navier-stokes :Re 100
                 '((lid :velocity (1 0))
                   (bottom :no-slip)
                   (left :no-slip)
                   (right :no-slip)))

               (sdl:time :dt 0.001 :end 10.0)

               (sdl:output :vtk :every 0.1))))
    (check-not-null sim "Simulation created")
    (check-equal (sdl:simulation-name sim) "lid-driven-cavity" "Simulation name")
    (check-equal (sdl:simulation-dim sim) 2 "Simulation dimension")
    (check-not-null (sdl:simulation-domain sim) "Domain defined")
    (check-not-null (sdl:simulation-boundaries sim) "Boundaries defined")
    (check-not-null (sdl:simulation-physics sim) "Physics defined")))

(define-test test-simulation-multiphysics
  "Test multiphysics simulation definition."
  (let ((sim (sdl:simulation "conjugate-heat" :dim 2
               (sdl:domain (sdl:box '(0 0) '(2 1)) :mesh :n '(80 40))

               (sdl:subdomains
                 '((fluid (< x 1))
                   (solid (>= x 1))))

               (sdl:boundaries
                 '((inlet (= x 0))
                   (outlet (= x 2))
                   (walls (or (= y 0) (= y 1)))))

               (sdl:physics :flow :navier-stokes
                 :subdomain 'fluid
                 :Re 100
                 '((inlet :velocity (0.1 0))
                   (outlet :pressure 0)
                   (walls :no-slip)))

               (sdl:physics :thermal :heat-conduction
                 '((inlet :temperature 20)
                   (outlet :outflow)
                   (walls :insulated)))

               (sdl:coupling :flow-thermal
                 :physics '(flow thermal)
                 :type :sequential))))
    (check-not-null sim "Multiphysics simulation created")
    (check-equal (length (sdl:simulation-subdomains sim)) 2 "Two subdomains")
    (check-equal (length (sdl:simulation-physics-list sim)) 2 "Two physics")
    (check-not-null (sdl:simulation-coupling sim) "Coupling defined")))

;;; ============================================================
;;; Run All Tests
;;; ============================================================

(defun run-all-tests ()
  "Run all SDL tests."
  (setf *test-count* 0
        *pass-count* 0
        *fail-count* 0)

  (format t "~%========================================~%")
  (format t "StreamVorti SDL Test Suite~%")
  (format t "========================================~%")

  ;; Geometry tests
  (format t "~%--- Geometry Primitives ---~%")
  (test-box-creation)
  (test-ball-creation)
  (test-geometry-contains)

  ;; Domain tests
  (format t "~%--- Domain Creation ---~%")
  (test-domain-mesh-structured)
  (test-domain-points)
  (test-domain-from-file)
  (test-domain-particles)

  ;; Boundary tests
  (format t "~%--- Boundaries ---~%")
  (test-boundary-predicate-equality)
  (test-boundary-predicate-boolean)
  (test-boundary-attribute)
  (test-boundaries-block)

  ;; Subdomain tests
  (format t "~%--- Subdomains ---~%")
  (test-subdomains-block)

  ;; Physics tests
  (format t "~%--- Physics ---~%")
  (test-physics-navier-stokes)
  (test-physics-heat-conduction)
  (test-boundary-conditions)
  (test-bc-with-function)

  ;; Method tests
  (format t "~%--- Methods ---~%")
  (test-method-dcpse)
  (test-method-fem)
  (test-method-sph)

  ;; Solver tests
  (format t "~%--- Solvers ---~%")
  (test-time-solver)
  (test-steady-solver)
  (test-linear-solver)

  ;; Material tests
  (format t "~%--- Materials ---~%")
  (test-material-definition)
  (test-fluid-material)

  ;; Output tests
  (format t "~%--- Output ---~%")
  (test-output-basic)
  (test-output-probes)

  ;; Coupling tests
  (format t "~%--- Coupling ---~%")
  (test-coupling-sequential)
  (test-coupling-monolithic)

  ;; Complete simulation tests
  (format t "~%--- Complete Simulations ---~%")
  (test-simulation-lid-driven-cavity)
  (test-simulation-multiphysics)

  ;; Summary
  (format t "~%========================================~%")
  (format t "Results: ~D passed, ~D failed, ~D total~%"
          *pass-count* *fail-count* *test-count*)
  (format t "========================================~%")

  (zerop *fail-count*))

;;; Run tests when loaded
(unless (run-all-tests)
  (format t "~%TESTS FAILED~%"))
