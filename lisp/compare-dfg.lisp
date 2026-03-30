;;;; compare-dfg.lisp - Validation against DFG cylinder benchmark
;;;;
;;;; Copyright (c) 2026 Benjamin F. Zwick
;;;;
;;;; Compares StreamVorti cylinder flow results against:
;;;;   Schäfer & Turek (1996), "Benchmark Computations of Laminar
;;;;   Flow Around a Cylinder", Flow Simulation with High-Performance
;;;;   Computers II, Notes on Numerical Fluid Mechanics, vol. 52.
;;;;
;;;; Usage:
;;;;   sbcl --load lisp/compare-dfg.lisp \
;;;;        --eval '(streamvorti.dfg:run-validation
;;;;                 :results-file "output_dat/dfg-cylinder_mid-channel.dat"
;;;;                 :H 0.41 :U-max 0.3)'

(defpackage :streamvorti.dfg
  (:use :cl)
  (:export #:run-validation))

(in-package :streamvorti.dfg)

;;; ============================================================
;;; Utilities (shared with compare-ghia.lisp pattern)
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
      (when (numberp val) (coerce val 'double-float)))))

(defun load-centerline (path)
  "Load centerline data file. Returns list of (y u v) triples."
  (unless (probe-file path)
    (error "Results file not found: ~A" path))
  (with-open-file (s path :direction :input)
    (loop for line = (read-line s nil nil)
          while line
          when (and (plusp (length line))
                    (char/= (char line 0) #\#))
            collect (let ((parts (split-whitespace line)))
                      (list (parse-number (first parts))
                            (parse-number (second parts))
                            (parse-number (third parts)))))))

(defun load-reference (path)
  "Load DFG reference data file. Returns alist of (name . value)."
  (unless (probe-file path)
    (error "Reference file not found: ~A" path))
  (with-open-file (s path :direction :input)
    (loop for line = (read-line s nil nil)
          while line
          when (and (plusp (length line))
                    (char/= (char line 0) #\#))
            collect (let ((parts (split-whitespace line)))
                      (cons (first parts)
                            (parse-number (second parts)))))))

(defun trapz (data key-x key-y)
  "Trapezoidal integration of KEY-Y over KEY-X from DATA points."
  (loop for i from 1 below (length data)
        sum (let ((x0 (funcall key-x (nth (1- i) data)))
                  (x1 (funcall key-x (nth i data)))
                  (y0 (funcall key-y (nth (1- i) data)))
                  (y1 (funcall key-y (nth i data))))
              (* 0.5d0 (+ y0 y1) (- x1 x0)))))

;;; ============================================================
;;; Validation
;;; ============================================================

(defun run-validation (&key results-file
                            (reference-file "data/turek_1996_dfg_re20.txt")
                            (H 0.41d0)
                            (U-max 0.3d0))
  "Compare StreamVorti results with DFG benchmark reference."
  (format t "~%======================================================================~%")
  (format t "DFG Cylinder Benchmark Validation~%")
  (format t "======================================================================~%")
  (format t "Reference: Schäfer & Turek (1996)~%~%")

  (unless results-file
    (setf results-file "output_dat/dfg-cylinder_mid-channel.dat"))

  ;; Load data
  (let ((data (load-centerline results-file))
        (ref (when (probe-file reference-file)
               (load-reference reference-file))))

    (format t "Loaded ~D data points from ~A~%" (length data) results-file)
    (when ref
      (format t "Reference values from ~A:~%" reference-file)
      (dolist (r ref)
        (format t "  ~A = ~,12F~%" (car r) (cdr r)))
      (format t "~%"))

    ;; Velocity profile table
    (format t "Velocity profile at probe location:~%")
    (format t "~8A  ~12A  ~12A~%" "y" "u" "v")
    (format t "~8A  ~12A  ~12A~%" "--------" "------------" "------------")
    (dolist (pt data)
      (format t "~8,4F  ~12,6F  ~12,6F~%" (first pt) (second pt) (third pt)))

    ;; Mass conservation: ∫u dy = Q = 2/3 · U_max · H
    (let* ((Q-ref (* 2/3 U-max H))
           (Q-computed (trapz data #'first #'second))
           (Q-error (if (> (abs Q-ref) 1d-15)
                        (* 100.0d0 (abs (/ (- Q-computed Q-ref) Q-ref)))
                        0.0d0)))
      (format t "~%Mass conservation:~%")
      (format t "  Inlet flow rate Q = 2/3·U_max·H = ~,6F~%" Q-ref)
      (format t "  ∫u dy at probe    = ~,6F~%" Q-computed)
      (format t "  Relative error    = ~,2F%~%" Q-error)
      (when (> Q-error 10.0)
        (format t "  WARNING: mass conservation error > 10%~%")))

    ;; Velocity extremes
    (let ((u-max (reduce #'max data :key #'second))
          (u-min (reduce #'min data :key #'second))
          (v-max (reduce #'max data :key #'third))
          (v-min (reduce #'min data :key #'third)))
      (format t "~%Velocity extremes:~%")
      (format t "  u: [~,6F, ~,6F]~%" u-min u-max)
      (format t "  v: [~,6F, ~,6F]~%" v-min v-max))

    ;; ASCII plot
    (format t "~%u-velocity profile:~%")
    (let ((u-max (max 1d-15 (reduce #'max data :key (lambda (p) (abs (second p)))))))
      (dolist (pt data)
        (let* ((y (first pt))
               (u (second pt))
               (bar (round (* 50 (/ (abs u) u-max))))
               (ch (if (>= u 0) #\# #\-)))
          (format t "  y=~5,3F |~A (~,4F)~%"
                  y (make-string (max 0 bar) :initial-element ch) u))))

    (format t "~%======================================================================~%")))
