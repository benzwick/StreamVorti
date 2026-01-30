;;;; compare-ghia.lisp - Compare StreamVorti results with Ghia et al. (1982)
;;;;
;;;; This script loads centerline data from StreamVorti and compares it
;;;; against the benchmark data from Ghia et al. (1982).
;;;;
;;;; Usage from SBCL/ECL REPL:
;;;;   (load "lisp/compare-ghia.lisp")
;;;;   (in-package :streamvorti.validation)
;;;;   (run-validation :results-file "output_dat/cavity_u_centerline_x0.5.dat")
;;;;
;;;; Or from command line:
;;;;   sbcl --load lisp/compare-ghia.lisp \
;;;;        --eval '(streamvorti.validation:run-validation)'
;;;;

(defpackage :streamvorti.validation
  (:use :cl)
  (:export #:load-centerline-data
           #:load-ghia-data
           #:compute-errors
           #:ascii-plot
           #:run-validation))

(in-package :streamvorti.validation)

;;; ============================================================
;;; Reference data for validation
;;; ============================================================

;; Ghia et al. (1982) - u-velocity along vertical centerline (x=0.5)
;; Format: (y . u-velocity) for Re=100
(defparameter *ghia-re100-u*
  '((1.0000 . 1.00000)
    (0.9766 . 0.84123)
    (0.9688 . 0.78871)
    (0.9609 . 0.73722)
    (0.9531 . 0.68717)
    (0.8516 . 0.23151)
    (0.7344 . 0.00332)
    (0.6172 . -0.13641)
    (0.5000 . -0.20581)
    (0.4531 . -0.21090)
    (0.2813 . -0.15662)
    (0.1719 . -0.10150)
    (0.1016 . -0.06434)
    (0.0703 . -0.04775)
    (0.0625 . -0.04192)
    (0.0547 . -0.03717)
    (0.0000 . 0.00000))
  "Ghia et al. (1982) u-velocity at x=0.5 for Re=100")

(defparameter *ghia-re400-u*
  '((1.0000 . 1.00000)
    (0.9766 . 0.75837)
    (0.9688 . 0.68439)
    (0.9609 . 0.61756)
    (0.9531 . 0.55892)
    (0.8516 . 0.29093)
    (0.7344 . 0.16257)
    (0.6172 . 0.02135)
    (0.5000 . -0.11477)
    (0.4531 . -0.17119)
    (0.2813 . -0.32726)
    (0.1719 . -0.24299)
    (0.1016 . -0.14612)
    (0.0703 . -0.10338)
    (0.0625 . -0.09266)
    (0.0547 . -0.08186)
    (0.0000 . 0.00000))
  "Ghia et al. (1982) u-velocity at x=0.5 for Re=400")

(defparameter *ghia-re1000-u*
  '((1.0000 . 1.00000)
    (0.9766 . 0.65928)
    (0.9688 . 0.57492)
    (0.9609 . 0.51117)
    (0.9531 . 0.46604)
    (0.8516 . 0.33304)
    (0.7344 . 0.18719)
    (0.6172 . 0.05702)
    (0.5000 . -0.06080)
    (0.4531 . -0.10648)
    (0.2813 . -0.27805)
    (0.1719 . -0.38289)
    (0.1016 . -0.29730)
    (0.0703 . -0.22220)
    (0.0625 . -0.20196)
    (0.0547 . -0.18109)
    (0.0000 . 0.00000))
  "Ghia et al. (1982) u-velocity at x=0.5 for Re=1000")

;;; ============================================================
;;; Data loading functions
;;; ============================================================

(defun parse-float-safe (string)
  "Parse a floating point number from STRING, returning NIL on error."
  (handler-case
      (let ((result (read-from-string string)))
        (if (numberp result) (float result) nil))
    (error () nil)))

(defun split-whitespace (string)
  "Split STRING by whitespace into a list of substrings."
  (let ((result nil)
        (start 0)
        (len (length string)))
    (loop
      (loop while (and (< start len)
                       (member (char string start) '(#\Space #\Tab)))
            do (incf start))
      (when (>= start len)
        (return (nreverse result)))
      (let ((end (position-if (lambda (c) (member c '(#\Space #\Tab #\Newline)))
                              string :start start)))
        (push (subseq string start (or end len)) result)
        (setf start (or end len))))))

(defun load-centerline-data (filename)
  "Load centerline data from file. Returns list of (y . u) pairs."
  (with-open-file (stream filename :direction :input :if-does-not-exist nil)
    (unless stream
      (format t "Error: Cannot open file ~A~%" filename)
      (return-from load-centerline-data nil))
    (loop for line = (read-line stream nil nil)
          while line
          unless (or (zerop (length line))
                     (char= (char line 0) #\#))
            collect (let* ((parts (split-whitespace line))
                           (y (parse-float-safe (first parts)))
                           (u (parse-float-safe (second parts))))
                      (when (and y u)
                        (cons y u))))))

(defun load-ghia-data (reynolds)
  "Get Ghia reference data for given Reynolds number."
  (case reynolds
    (100 *ghia-re100-u*)
    (400 *ghia-re400-u*)
    (1000 *ghia-re1000-u*)
    (otherwise
     (format t "Warning: No Ghia data for Re=~A, using Re=100~%" reynolds)
     *ghia-re100-u*)))

;;; ============================================================
;;; Error computation functions
;;; ============================================================

(defun interpolate (x-target data)
  "Linear interpolation of DATA (sorted by car) at X-TARGET."
  (let ((sorted (sort (copy-list data) #'< :key #'car)))
    (cond
      ((<= x-target (caar sorted)) (cdar sorted))
      ((>= x-target (caar (last sorted))) (cdar (last sorted)))
      (t (loop for (p1 p2) on sorted
               when (and p2 (<= (car p1) x-target (car p2)))
                 return (let ((t-val (/ (- x-target (car p1))
                                        (- (car p2) (car p1)))))
                          (+ (cdr p1) (* t-val (- (cdr p2) (cdr p1))))))))))

(defun compute-errors (computed-data reference-data)
  "Compute L2 and Linf errors between computed and reference data.
   Returns (values l2-error linf-error)."
  (let* ((errors (mapcar (lambda (ref-point)
                           (let ((y (car ref-point))
                                 (u-ref (cdr ref-point)))
                             (- (interpolate y computed-data) u-ref)))
                         reference-data))
         (squared-errors (mapcar (lambda (e) (* e e)) errors))
         (l2-error (sqrt (/ (reduce #'+ squared-errors) (length errors))))
         (linf-error (reduce #'max (mapcar #'abs errors))))
    (values l2-error linf-error)))

;;; ============================================================
;;; ASCII plotting (no dependencies)
;;; ============================================================

(defun ascii-plot (computed-data reference-data &key (width 60) (height 20))
  "Simple ASCII plot for terminal output.
   Plots u-velocity (x-axis) vs y (y-axis)."
  ;; Filter out nil entries
  (setf computed-data (remove-if #'null computed-data))
  (setf reference-data (remove-if #'null reference-data))

  (let* ((all-u (append (mapcar #'cdr computed-data) (mapcar #'cdr reference-data)))
         (u-min (reduce #'min all-u))
         (u-max (reduce #'max all-u))
         (u-range (- u-max u-min))
         (canvas (make-array (list height width) :initial-element #\Space)))

    ;; Plot reference points as 'o'
    (dolist (pt reference-data)
      (when pt
        (let ((row (floor (* (- 1.0 (car pt)) (1- height))))
              (col (floor (* (/ (- (cdr pt) u-min) u-range) (1- width)))))
          (when (and (<= 0 row (1- height)) (<= 0 col (1- width)))
            (setf (aref canvas row col) #\o)))))

    ;; Plot computed points as '*'
    (dolist (pt computed-data)
      (when pt
        (let ((row (floor (* (- 1.0 (car pt)) (1- height))))
              (col (floor (* (/ (- (cdr pt) u-min) u-range) (1- width)))))
          (when (and (<= 0 row (1- height)) (<= 0 col (1- width)))
            (setf (aref canvas row col) #\*)))))

    ;; Print header
    (format t "~%u-velocity vs y (o=Ghia, *=computed)~%")
    (format t "u: ~,2f~vT~,2f~%" u-min width u-max)

    ;; Print canvas
    (dotimes (row height)
      (format t "|")
      (dotimes (col width)
        (format t "~c" (aref canvas row col)))
      (format t "|~%"))

    (format t "y: 1.0~vT0.0~%" (- width 3))))

;;; ============================================================
;;; Main validation function
;;; ============================================================

(defun run-validation (&key
                         (results-file "output_dat/mfem_square10x10_u_centerline_x0.5.dat")
                         (reynolds 100))
  "Run full validation: load data, compute errors, generate ASCII plot.

   Keyword arguments:
     :results-file  Path to StreamVorti centerline data file
     :reynolds      Reynolds number (100, 400, or 1000 for Ghia data)"

  (format t "~%~A~%" (make-string 70 :initial-element #\=))
  (format t "StreamVorti Validation vs Ghia et al. (1982)~%")
  (format t "~A~%~%" (make-string 70 :initial-element #\=))

  (format t "Reynolds number: ~A~%" reynolds)
  (format t "Results file: ~A~%" results-file)

  ;; Load computed results
  (let ((computed (load-centerline-data results-file)))
    (unless computed
      (format t "~%ERROR: Could not load results file~%")
      (return-from run-validation nil))

    (format t "~%Loaded ~A computed data points~%" (length computed))

    ;; Get reference data
    (let ((reference (load-ghia-data reynolds)))
      (format t "Reference data: ~A points~%" (length reference))

      ;; Compute errors
      (multiple-value-bind (l2 linf) (compute-errors computed reference)
        (format t "~%Error metrics:~%")
        (format t "  L2 error:   ~,6f~%" l2)
        (format t "  Linf error: ~,6f~%" linf)

        ;; Pass/fail criteria
        (let ((l2-threshold 0.05)
              (linf-threshold 0.10))
          (format t "~%Validation result:~%")
          (format t "  L2 error < ~,0f%: ~A~%"
                  (* 100 l2-threshold) (if (< l2 l2-threshold) "PASS" "FAIL"))
          (format t "  Linf error < ~,0f%: ~A~%"
                  (* 100 linf-threshold) (if (< linf linf-threshold) "PASS" "FAIL"))))

      ;; Generate ASCII plot
      (format t "~%")
      (ascii-plot computed reference)

      (format t "~%~A~%" (make-string 70 :initial-element #\=))

      ;; Return the computed data for further analysis
      computed)))

;;; Print usage when loaded
(format t "~%StreamVorti Validation Script Loaded~%")
(format t "Usage: (streamvorti.validation:run-validation :results-file \"path/to/centerline.dat\")~%~%")
