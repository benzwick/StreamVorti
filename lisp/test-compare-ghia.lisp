;;;; test-compare-ghia.lisp - Test suite for compare-ghia.lisp
;;;;
;;;; Run with SBCL:
;;;;   sbcl --load lisp/compare-ghia.lisp --load lisp/test-compare-ghia.lisp
;;;;
;;;; Run with ECL:
;;;;   ecl -load lisp/compare-ghia.lisp -load lisp/test-compare-ghia.lisp

(defpackage :streamvorti.validation.tests
  (:use :cl :streamvorti.validation)
  (:export #:run-all-tests))

(in-package :streamvorti.validation.tests)

;;; ============================================================
;;; Test Framework (minimal, no dependencies)
;;; ============================================================

(defvar *test-count* 0)
(defvar *pass-count* 0)
(defvar *fail-count* 0)
(defvar *test-results* nil)

(defmacro define-test (name &body body)
  "Define a test case."
  `(defun ,name ()
     (handler-case
         (progn ,@body t)
       (error (e)
         (format t "~%ERROR in ~A: ~A~%" ',name e)
         nil))))

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

;;; ============================================================
;;; Tests for parsing utilities
;;; ============================================================

(define-test test-parse-float-safe
  "Test parse-float-safe with various inputs."
  (check-near (streamvorti.validation::parse-float-safe "3.14") 3.14 0.001
              "Parse positive float")
  (check-near (streamvorti.validation::parse-float-safe "-0.5") -0.5 0.001
              "Parse negative float")
  (check-near (streamvorti.validation::parse-float-safe "42") 42.0 0.001
              "Parse integer as float")
  (check-near (streamvorti.validation::parse-float-safe "1.0e-3") 0.001 1e-6
              "Parse scientific notation")
  (check (null (streamvorti.validation::parse-float-safe "abc"))
         "Return NIL for non-numeric string")
  (check (null (streamvorti.validation::parse-float-safe ""))
         "Return NIL for empty string"))

(define-test test-split-whitespace
  "Test split-whitespace with various inputs."
  (check-equal (streamvorti.validation::split-whitespace "a b c")
               '("a" "b" "c")
               "Split simple space-separated")
  (check-equal (streamvorti.validation::split-whitespace "  a   b  c  ")
               '("a" "b" "c")
               "Handle multiple spaces")
  (check-equal (streamvorti.validation::split-whitespace "a	b	c")
               '("a" "b" "c")
               "Handle tabs")
  (check-equal (streamvorti.validation::split-whitespace "")
               nil
               "Empty string returns nil")
  (check-equal (streamvorti.validation::split-whitespace "   ")
               nil
               "Whitespace-only returns nil")
  (check-equal (streamvorti.validation::split-whitespace "single")
               '("single")
               "Single word"))

(define-test test-parse-reynolds-from-header
  "Test parsing Reynolds numbers from header line."
  (check-equal (streamvorti.validation::parse-reynolds-from-header
                "# y Re=100 Re=400 Re=1000")
               '(100 400 1000)
               "Parse standard header")
  (check-equal (streamvorti.validation::parse-reynolds-from-header
                "# y      Re=100   Re=400   Re=1000  Re=3200")
               '(100 400 1000 3200)
               "Parse with extra spaces")
  (check-equal (streamvorti.validation::parse-reynolds-from-header
                "# no reynolds numbers here")
               nil
               "Return nil when no Re= patterns")
  (check-equal (streamvorti.validation::parse-reynolds-from-header
                "Re=100")
               '(100)
               "Parse single Re without prefix"))

;;; ============================================================
;;; Tests for interpolation
;;; ============================================================

(define-test test-interpolate
  "Test linear interpolation."
  (let ((data '((0.0 . 0.0) (1.0 . 1.0))))
    (check-near (streamvorti.validation::interpolate 0.5 data) 0.5 1e-6
                "Interpolate midpoint")
    (check-near (streamvorti.validation::interpolate 0.0 data) 0.0 1e-6
                "Interpolate at start")
    (check-near (streamvorti.validation::interpolate 1.0 data) 1.0 1e-6
                "Interpolate at end")
    (check-near (streamvorti.validation::interpolate 0.25 data) 0.25 1e-6
                "Interpolate quarter point"))

  (let ((data '((0.0 . 10.0) (1.0 . 20.0) (2.0 . 30.0))))
    (check-near (streamvorti.validation::interpolate 0.5 data) 15.0 1e-6
                "Interpolate three points - first segment")
    (check-near (streamvorti.validation::interpolate 1.5 data) 25.0 1e-6
                "Interpolate three points - second segment"))

  ;; Test extrapolation (clamps to boundary values)
  (let ((data '((0.0 . 5.0) (1.0 . 10.0))))
    (check-near (streamvorti.validation::interpolate -1.0 data) 5.0 1e-6
                "Extrapolate below range returns first value")
    (check-near (streamvorti.validation::interpolate 2.0 data) 10.0 1e-6
                "Extrapolate above range returns last value")))

;;; ============================================================
;;; Tests for error computation
;;; ============================================================

(define-test test-compute-errors
  "Test L2 and Linf error computation."
  ;; Perfect match
  (let ((computed '((0.0 . 1.0) (0.5 . 2.0) (1.0 . 3.0)))
        (reference '((0.0 . 1.0) (0.5 . 2.0) (1.0 . 3.0))))
    (multiple-value-bind (l2 linf) (compute-errors computed reference)
      (check-near l2 0.0 1e-6 "L2 error for identical data")
      (check-near linf 0.0 1e-6 "Linf error for identical data")))

  ;; Constant offset
  (let ((computed '((0.0 . 1.1) (0.5 . 2.1) (1.0 . 3.1)))
        (reference '((0.0 . 1.0) (0.5 . 2.0) (1.0 . 3.0))))
    (multiple-value-bind (l2 linf) (compute-errors computed reference)
      (check-near l2 0.1 1e-6 "L2 error for constant offset")
      (check-near linf 0.1 1e-6 "Linf error for constant offset")))

  ;; Single large error
  (let ((computed '((0.0 . 1.0) (0.5 . 2.5) (1.0 . 3.0)))
        (reference '((0.0 . 1.0) (0.5 . 2.0) (1.0 . 3.0))))
    (multiple-value-bind (l2 linf) (compute-errors computed reference)
      (check-near linf 0.5 1e-6 "Linf error captures max error")
      (check (< l2 linf) "L2 error less than Linf for single outlier"))))

;;; ============================================================
;;; Tests for file loading
;;; ============================================================

(define-test test-load-reference-file
  "Test loading reference data files."
  ;; Test loading Ghia data
  (multiple-value-bind (reynolds data-hash)
      (streamvorti.validation::load-reference-file "ghia_1982_u.txt")
    (check reynolds "Reynolds list should not be nil")
    (check data-hash "Data hash should not be nil")
    (when reynolds
      (check (member 100 reynolds) "Should have Re=100")
      (check (member 1000 reynolds) "Should have Re=1000")
      (check (member 10000 reynolds) "Should have Re=10000"))
    (when data-hash
      (let ((re100-data (gethash 100 data-hash)))
        (check re100-data "Should have data for Re=100")
        (when re100-data
          (check (> (length re100-data) 10)
                 "Should have >10 data points for Re=100")
          ;; Check boundary values
          (let ((first-point (first re100-data))
                (last-point (car (last re100-data))))
            (check-near (car first-point) 1.0 0.01
                        "First y value should be 1.0")
            (check-near (cdr first-point) 1.0 0.01
                        "First u value should be 1.0 (lid velocity)")
            (check-near (car last-point) 0.0 0.01
                        "Last y value should be 0.0")
            (check-near (cdr last-point) 0.0 0.01
                        "Last u value should be 0.0 (no-slip)")))))))

(define-test test-load-reference-file-erturk
  "Test loading Erturk reference data."
  (multiple-value-bind (reynolds data-hash)
      (streamvorti.validation::load-reference-file "erturk_2005_u.txt")
    (declare (ignore data-hash))
    (check reynolds "Reynolds list should not be nil")
    (when reynolds
      (check (member 1000 reynolds) "Should have Re=1000")
      (check (member 5000 reynolds) "Should have Re=5000")
      (check (member 20000 reynolds) "Should have Re=20000"))))

(define-test test-load-reference-file-missing
  "Test handling of missing file."
  (multiple-value-bind (reynolds data-hash)
      (streamvorti.validation::load-reference-file "nonexistent_file.txt")
    (check (null reynolds) "Reynolds should be nil for missing file")
    (check (null data-hash) "Data hash should be nil for missing file")))

(define-test test-load-centerline-data
  "Test loading centerline data from simulation output."
  ;; Create a temporary test file
  (let ((test-file "/tmp/test_centerline.dat"))
    (with-open-file (out test-file :direction :output :if-exists :supersede)
      (format out "# Test centerline data~%")
      (format out "# y  u  v~%")
      (format out "0.0 0.0 0.0~%")
      (format out "0.5 -0.1 0.05~%")
      (format out "1.0 1.0 0.0~%"))
    (let ((data (load-centerline-data test-file)))
      (check data "Should load data")
      (when data
        (check-equal (length data) 3 "Should have 3 data points")
        (check-near (car (first data)) 0.0 1e-6 "First y = 0.0")
        (check-near (cdr (first data)) 0.0 1e-6 "First u = 0.0")
        (check-near (car (second data)) 0.5 1e-6 "Second y = 0.5")
        (check-near (cdr (second data)) -0.1 1e-6 "Second u = -0.1")
        (check-near (car (third data)) 1.0 1e-6 "Third y = 1.0")
        (check-near (cdr (third data)) 1.0 1e-6 "Third u = 1.0")))
    ;; Cleanup
    (delete-file test-file)))

(define-test test-load-centerline-data-missing
  "Test handling of missing centerline file."
  (let ((data (load-centerline-data "/tmp/nonexistent_centerline.dat")))
    (check (null data) "Should return nil for missing file")))

;;; ============================================================
;;; Tests for get-reference-data
;;; ============================================================

(define-test test-load-reference-data
  "Test the high-level load-reference-data function."
  ;; Test Re=100 (Ghia only)
  (multiple-value-bind (data-list sources) (load-reference-data 100)
    (check data-list "Should have data for Re=100")
    (check sources "Should have sources for Re=100")
    (when sources
      (check (search "Ghia" (first sources))
             "Re=100 should come from Ghia")))

  ;; Test Re=1000 (both Ghia and Erturk)
  (multiple-value-bind (data-list sources) (load-reference-data 1000)
    (declare (ignore sources))
    (check data-list "Should have data for Re=1000")
    (check (>= (length data-list) 1) "Should have at least one source for Re=1000"))

  ;; Test Re=5000 (Erturk only for high Re)
  (multiple-value-bind (data-list sources) (load-reference-data 5000)
    (declare (ignore sources))
    (check data-list "Should have data for Re=5000")))

;;; ============================================================
;;; Tests for ASCII plot
;;; ============================================================

(define-test test-ascii-plot
  "Test ASCII plot generation (should not error)."
  (let ((computed '((0.0 . 0.0) (0.5 . -0.1) (1.0 . 1.0)))
        (reference '((0.0 . 0.0) (0.5 . -0.15) (1.0 . 1.0))))
    ;; Capture output to verify no errors
    (let ((output (with-output-to-string (*standard-output*)
                    (ascii-plot computed reference :width 40 :height 10))))
      (check (> (length output) 0) "ASCII plot should produce output")
      (check (search "u-velocity" output) "Output should contain axis label")
      (check (search "|" output) "Output should contain border characters"))))

(define-test test-ascii-plot-empty-data
  "Test ASCII plot handles empty/nil data gracefully."
  ;; Should not error with nil elements
  (let ((computed '((0.0 . 0.0) nil (1.0 . 1.0)))
        (reference '((0.0 . 0.0) (1.0 . 1.0))))
    (let ((output (with-output-to-string (*standard-output*)
                    (ascii-plot computed reference :width 30 :height 8))))
      (check (> (length output) 0) "Should handle nil elements"))))

;;; ============================================================
;;; Integration tests
;;; ============================================================

(define-test test-run-validation-with-test-data
  "Integration test with synthetic test data."
  ;; Create test data that matches Ghia Re=100 approximately
  (let ((test-file "/tmp/test_validation_centerline.dat"))
    (with-open-file (out test-file :direction :output :if-exists :supersede)
      (format out "# Test centerline~%")
      (format out "# y u v~%")
      ;; Approximate Ghia Re=100 values
      (format out "0.0000 0.0 0.0~%")
      (format out "0.0547 -0.037 0.0~%")
      (format out "0.0625 -0.042 0.0~%")
      (format out "0.0703 -0.048 0.0~%")
      (format out "0.1016 -0.064 0.0~%")
      (format out "0.1719 -0.102 0.0~%")
      (format out "0.2813 -0.157 0.0~%")
      (format out "0.4531 -0.211 0.0~%")
      (format out "0.5000 -0.206 0.0~%")
      (format out "0.6172 -0.136 0.0~%")
      (format out "0.7344 0.003 0.0~%")
      (format out "0.8516 0.232 0.0~%")
      (format out "0.9531 0.687 0.0~%")
      (format out "0.9609 0.737 0.0~%")
      (format out "0.9688 0.789 0.0~%")
      (format out "0.9766 0.841 0.0~%")
      (format out "1.0000 1.0 0.0~%"))

    ;; Run validation (capture output)
    (let ((output (with-output-to-string (*standard-output*)
                    (run-validation :results-file test-file :reynolds 100))))
      (check (search "Re=100" output) "Should mention Reynolds number")
      (check (search "Ghia" output) "Should mention Ghia reference")
      (check (search "L2 error" output) "Should report L2 error")
      (check (search "Linf error" output) "Should report Linf error"))

    ;; Cleanup
    (delete-file test-file)))

;;; ============================================================
;;; Test runner
;;; ============================================================

(defun run-all-tests ()
  "Run all tests and report results."
  (setf *test-count* 0
        *pass-count* 0
        *fail-count* 0)

  (format t "~%~A~%" (make-string 70 :initial-element #\=))
  (format t "StreamVorti Validation Test Suite~%")
  (format t "~A~%~%" (make-string 70 :initial-element #\=))

  (let ((tests '(test-parse-float-safe
                 test-split-whitespace
                 test-parse-reynolds-from-header
                 test-interpolate
                 test-compute-errors
                 test-load-reference-file
                 test-load-reference-file-erturk
                 test-load-reference-file-missing
                 test-load-centerline-data
                 test-load-centerline-data-missing
                 test-load-reference-data
                 test-ascii-plot
                 test-ascii-plot-empty-data
                 test-run-validation-with-test-data)))

    (dolist (test tests)
      (format t "~%Running ~A...~%" test)
      (let ((result (funcall test)))
        (if result
            (format t "  OK~%")
            (format t "  FAILED~%")))))

  (format t "~%~A~%" (make-string 70 :initial-element #\=))
  (format t "Test Results: ~A passed, ~A failed, ~A total checks~%"
          *pass-count* *fail-count* *test-count*)
  (format t "~A~%~%" (make-string 70 :initial-element #\=))

  (if (zerop *fail-count*)
      (format t "ALL TESTS PASSED~%")
      (format t "SOME TESTS FAILED~%"))

  (values *pass-count* *fail-count* *test-count*))

;;; Run tests when loaded
(run-all-tests)
