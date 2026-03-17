;;;; compare-poiseuille.lisp - Validation against exact analytical solution
;;;;
;;;; Copyright (c) 2026 Benjamin F. Zwick
;;;;
;;;; Compares StreamVorti Poiseuille/channel flow results against the
;;;; exact analytical solution for fully-developed flow between parallel
;;;; plates:
;;;;
;;;;   u(y) = 4 * U_max * y * (H - y) / H²
;;;;   v(y) = 0
;;;;
;;;; where H is the channel height and U_max is the maximum velocity
;;;; (at the centerline y = H/2).
;;;;
;;;; Usage:
;;;;   sbcl --load lisp/compare-poiseuille.lisp \
;;;;        --eval '(streamvorti.poiseuille:run-validation
;;;;                  :results-file "output_dat/poiseuille-flow_u_centerline_x0.5.dat"
;;;;                  :u-max 1.0 :channel-height 1.0)'

(defpackage :streamvorti.poiseuille
  (:use :cl)
  (:export #:run-validation
           #:load-centerline-data
           #:analytical-u
           #:compute-errors
           #:ascii-plot))

(in-package :streamvorti.poiseuille)

;;; ============================================================
;;; Utilities
;;; ============================================================

(defun split-whitespace (string)
  "Split STRING on whitespace."
  (loop with start = 0
        with len = (length string)
        for pos = (position-if (lambda (c) (member c '(#\Space #\Tab)))
                               string :start start)
        while (< start len)
        when (or (null pos) (> pos start))
          collect (subseq string start (or pos len))
        do (setf start (if pos (1+ pos) len))))

(defun parse-number (string)
  "Parse STRING as a number, returning NIL on failure."
  (ignore-errors
    (let ((val (read-from-string string nil nil)))
      (when (numberp val)
        (float val 1.0d0)))))

;;; ============================================================
;;; Analytical Solution
;;; ============================================================

(defun analytical-u (y &key (u-max 1.0d0) (channel-height 1.0d0))
  "Exact analytical u-velocity for Poiseuille flow.
   u(y) = 4 * U_max * y * (H - y) / H²"
  (let ((h channel-height))
    (* 4.0d0 u-max y (- h y) (/ 1.0d0 (* h h)))))

;;; ============================================================
;;; Simulation Results
;;; ============================================================

(defun load-centerline-data (filename)
  "Load centerline data from FILENAME. Returns list of (y . u) pairs."
  (unless (probe-file filename)
    (error "Results file not found: ~A" filename))
  (with-open-file (in filename)
    (loop for line = (read-line in nil)
          while line
          for tokens = (split-whitespace line)
          for y = (parse-number (first tokens))
          for u = (parse-number (second tokens))
          when (and y u (or (zerop (length line)) (char/= #\# (char line 0))))
            collect (cons y u))))

;;; ============================================================
;;; Error Computation
;;; ============================================================

(defun compute-errors (computed &key (u-max 1.0d0) (channel-height 1.0d0))
  "Compute L2 and L-infinity errors against exact Poiseuille solution.
Returns (VALUES l2 linf)."
  (let* ((errors (loop for (y . u-computed) in computed
                       for u-exact = (analytical-u y :u-max u-max
                                                     :channel-height channel-height)
                       collect (- u-computed u-exact)))
         (l2 (sqrt (/ (loop for e in errors sum (* e e)) (length errors))))
         (linf (loop for e in errors maximize (abs e))))
    (values l2 linf)))

;;; ============================================================
;;; ASCII Plot
;;; ============================================================

(defun ascii-plot (computed &key (u-max 1.0d0) (channel-height 1.0d0)
                              (width 60) (height 20))
  "Print ASCII plot comparing computed vs exact Poiseuille solution.
   Legend: o=exact, *=computed, @=both (overlapping)"
  ;; Generate exact solution at same y-values
  (let* ((exact (loop for (y . nil) in computed
                      collect (cons y (analytical-u y :u-max u-max
                                                      :channel-height channel-height))))
         (all-u (append (mapcar #'cdr computed) (mapcar #'cdr exact)))
         (u-min (reduce #'min all-u))
         (u-max-val (reduce #'max all-u))
         (u-range (max 1e-6 (- u-max-val u-min)))
         (canvas (make-array (list height width) :initial-element #\Space)))
    ;; Exact as 'o'
    (loop for (y . u) in exact
          for row = (min (1- height) (max 0 (floor (* (- 1.0 (/ y channel-height)) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- u u-min) u-range) (1- width)))))
          do (setf (aref canvas row col) #\o))
    ;; Computed as '*', or '@' if overlapping with exact
    (loop for (y . u) in computed
          for row = (min (1- height) (max 0 (floor (* (- 1.0 (/ y channel-height)) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- u u-min) u-range) (1- width)))))
          do (setf (aref canvas row col)
                   (if (char= (aref canvas row col) #\o) #\@ #\*)))
    ;; Print
    (format t "~%u-velocity vs y (o=exact, *=computed, @=both)~%")
    (format t "u: ~,2F~VT~,2F~%" u-min width u-max-val)
    (dotimes (row height)
      (format t "|")
      (dotimes (col width)
        (write-char (aref canvas row col)))
      (format t "|~%"))
    (format t "y: ~,1F~VT0.0~%" channel-height (- width 3))))

;;; ============================================================
;;; Main Entry Point
;;; ============================================================

(defun run-validation (&key
                         (results-file "output_dat/poiseuille-flow_u_centerline_x0.5.dat")
                         (u-max 1.0d0)
                         (channel-height 1.0d0)
                         (l2-threshold 0.05d0)
                         (linf-threshold 0.10d0))
  "Validate StreamVorti Poiseuille results against exact analytical solution.
Returns T if validation passes, NIL otherwise.

Exact solution: u(y) = 4 * U_max * y * (H - y) / H²

Keyword arguments:
  :results-file    - path to centerline data file
  :u-max           - maximum inlet velocity (default 1.0)
  :channel-height  - channel height H (default 1.0)
  :l2-threshold    - L2 error pass threshold (default 0.05 = 5%)
  :linf-threshold  - Linf error pass threshold (default 0.10 = 10%)"
  (format t "~%~70,,,'-<~>~%")
  (format t "StreamVorti Poiseuille Validation~%")
  (format t "~70,,,'-<~>~%~%")
  (format t "Results file: ~A~%" results-file)
  (format t "Parameters: U_max = ~A, H = ~A~%" u-max channel-height)
  (format t "Exact solution: u(y) = 4 * ~A * y * (~A - y) / ~A~%"
          u-max channel-height (* channel-height channel-height))
  ;; Load data
  (let ((computed (load-centerline-data results-file)))
    (format t "~%Loaded ~A computed data points~%" (length computed))
    ;; Compute errors
    (multiple-value-bind (l2 linf) (compute-errors computed :u-max u-max
                                                            :channel-height channel-height)
      (let ((l2-pass (< l2 l2-threshold))
            (linf-pass (< linf linf-threshold)))
        (format t "~%--- Exact Analytical Solution ---~%")
        (format t "  L2 error:   ~,6F~%" l2)
        (format t "  Linf error: ~,6F~%" linf)
        (format t "  L2 < ~A%: ~A~%" (round (* l2-threshold 100)) (if l2-pass "PASS" "FAIL"))
        (format t "  Linf < ~A%: ~A~%" (round (* linf-threshold 100)) (if linf-pass "PASS" "FAIL"))
        (format t "~%Overall: ~A~%" (if (and l2-pass linf-pass) "PASS" "FAIL"))
        ;; ASCII plot
        (ascii-plot computed :u-max u-max :channel-height channel-height)
        (format t "~%~70,,,'-<~>~%")
        (and l2-pass linf-pass)))))

;;; ============================================================

(format t "~%StreamVorti Poiseuille Validation Script Loaded~%")
(format t "Usage: (streamvorti.poiseuille:run-validation :results-file \"path/to/centerline.dat\")~%~%")
