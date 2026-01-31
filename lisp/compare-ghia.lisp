;;;; compare-ghia.lisp - Compare StreamVorti results with reference data
;;;;
;;;; Compares simulation results against:
;;;; - Ghia, Ghia & Shin (1982) for Re <= 1000
;;;; - Erturk & Corke (2005) for Re > 1000
;;;;
;;;; Reference data is loaded from data/*.txt files (deal.II format).
;;;;
;;;; Usage from SBCL/ECL REPL:
;;;;   (load "lisp/compare-ghia.lisp")
;;;;   (in-package :streamvorti.validation)
;;;;   (run-validation :results-file "output_dat/cavity_u_centerline_x0.5.dat")
;;;;
;;;; Or from command line:
;;;;   sbcl --load lisp/compare-ghia.lisp \
;;;;        --eval '(streamvorti.validation:run-validation)'

(defpackage :streamvorti.validation
  (:use :cl)
  (:export #:load-centerline-data
           #:load-reference-data
           #:compute-errors
           #:generate-plot
           #:ascii-plot
           #:run-validation))

(in-package :streamvorti.validation)

;;; ============================================================
;;; Configuration
;;; ============================================================

(defparameter *data-dir* "data/"
  "Directory containing reference data files.")

;;; ============================================================
;;; Parsing utilities
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

;;; ============================================================
;;; Loading reference data from files
;;; ============================================================

(defun parse-reynolds-from-header (header-line)
  "Parse Reynolds numbers from header line like '# y Re=100 Re=400 Re=1000'.
   Returns list of Reynolds numbers as integers."
  (let ((parts (split-whitespace header-line))
        (reynolds-list nil))
    (dolist (part parts)
      (when (and (> (length part) 3)
                 (string-equal "Re=" (subseq part 0 3)))
        (let ((re (parse-integer (subseq part 3) :junk-allowed t)))
          (when re (push re reynolds-list)))))
    (nreverse reynolds-list)))

(defun load-reference-file (filename)
  "Load multi-column reference data file.
   Returns (values reynolds-list data-hash) where data-hash maps Re -> ((y . u) ...)."
  (let ((filepath (merge-pathnames filename *data-dir*))
        (reynolds-list nil)
        (data-hash (make-hash-table)))
    (with-open-file (stream filepath :direction :input :if-does-not-exist nil)
      (unless stream
        (format t "Warning: Cannot open ~A~%" filepath)
        (return-from load-reference-file (values nil nil)))
      ;; Read header lines to get Reynolds numbers
      (loop for line = (read-line stream nil nil)
            while (and line (> (length line) 0) (char= (char line 0) #\#))
            do (let ((re-list (parse-reynolds-from-header line)))
                 (when re-list (setf reynolds-list re-list))))
      ;; Initialize hash table entries
      (dolist (re reynolds-list)
        (setf (gethash re data-hash) nil))
      ;; Rewind and read data
      (file-position stream 0)
      (loop for line = (read-line stream nil nil)
            while line
            unless (or (zerop (length line)) (char= (char line 0) #\#))
              do (let* ((parts (split-whitespace line))
                        (y (parse-float-safe (first parts))))
                   (when y
                     (loop for re in reynolds-list
                           for u-str in (rest parts)
                           for u = (parse-float-safe u-str)
                           when u do (push (cons y u) (gethash re data-hash)))))))
    ;; Reverse the lists (they were built in reverse order)
    (dolist (re reynolds-list)
      (setf (gethash re data-hash) (nreverse (gethash re data-hash))))
    (values reynolds-list data-hash)))

;;; Cache loaded data
(defvar *ghia-data* nil)
(defvar *ghia-reynolds* nil)
(defvar *erturk-data* nil)
(defvar *erturk-reynolds* nil)

(defun ensure-reference-data-loaded ()
  "Load reference data files if not already loaded."
  (unless *ghia-data*
    (multiple-value-setq (*ghia-reynolds* *ghia-data*)
      (load-reference-file "ghia_1982_u.txt")))
  (unless *erturk-data*
    (multiple-value-setq (*erturk-reynolds* *erturk-data*)
      (load-reference-file "erturk_2005_u.txt"))))

(defun load-reference-data (reynolds)
  "Get reference data for given Reynolds number.
   Returns data from both Ghia (1982) and Erturk (2005) when available.
   Returns (values data-list source-names) where data-list and source-names are lists."
  (ensure-reference-data-loaded)
  (let ((ghia-data (and *ghia-data* (gethash reynolds *ghia-data*)))
        (erturk-data (and *erturk-data* (gethash reynolds *erturk-data*)))
        (data-list nil)
        (source-names nil))
    ;; Collect available data from both sources
    (when ghia-data
      (push ghia-data data-list)
      (push "Ghia et al. (1982)" source-names))
    (when erturk-data
      (push erturk-data data-list)
      (push "Erturk & Corke (2005)" source-names))
    ;; If no data found, fall back to closest available
    (unless data-list
      (format t "Warning: No reference data for Re=~A~%" reynolds)
      (cond
        ((<= reynolds 1000)
         (push (gethash 100 *ghia-data*) data-list)
         (push "Ghia et al. (1982) [Re=100]" source-names))
        (t
         (push (gethash 1000 *erturk-data*) data-list)
         (push "Erturk & Corke (2005) [Re=1000]" source-names))))
    (values (nreverse data-list) (nreverse source-names))))

;;; ============================================================
;;; Loading simulation results
;;; ============================================================

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

;;; ============================================================
;;; Error computation
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
;;; Plotting with vgplot
;;; ============================================================

(defvar *vgplot-available* nil
  "Flag indicating if vgplot is loaded and available.")

(defun ensure-vgplot ()
  "Load vgplot via quicklisp if available. Returns T if successful."
  (when *vgplot-available*
    (return-from ensure-vgplot t))
  (handler-case
      (progn
        ;; Load quicklisp if not already loaded
        (unless (find-package :ql)
          (let ((quicklisp-init (merge-pathnames "quicklisp/setup.lisp"
                                                  (user-homedir-pathname))))
            (when (probe-file quicklisp-init)
              (load quicklisp-init))))
        ;; Load vgplot via quicklisp
        (when (find-package :ql)
          (funcall (intern "QUICKLOAD" :ql) :vgplot :silent t)
          (setf *vgplot-available* t)
          t))
    (error (e)
      (format t "Note: vgplot not available (~A)~%" e)
      nil)))

(defun generate-plot (computed-data reference-data reynolds source-name
                      &key (output-dir "plots/") (output-file "validation_plot.png"))
  "Generate publication-quality plot using vgplot (Common Lisp gnuplot interface).
   Returns path to generated PNG file, or NIL if vgplot not available."
  (unless (ensure-vgplot)
    (return-from generate-plot nil))

  (let ((vgplot (find-package :vgplot)))
    (unless vgplot
      (return-from generate-plot nil))

    ;; Ensure output directory exists
    (ensure-directories-exist (make-pathname :directory (pathname-directory
                                                          (merge-pathnames output-file output-dir))))

    (let ((png-path (merge-pathnames output-file output-dir)))
      ;; Extract data as lists
      (let ((ref-y (mapcar #'car (remove-if #'null reference-data)))
            (ref-u (mapcar #'cdr (remove-if #'null reference-data)))
            (comp-y (mapcar #'car (remove-if #'null computed-data)))
            (comp-u (mapcar #'cdr (remove-if #'null computed-data)))
            (format-plot (intern "FORMAT-PLOT" vgplot)))

        ;; Set plot parameters
        (funcall (intern "TITLE" vgplot)
                 (format nil "Lid-Driven Cavity: u-velocity at x=0.5 (Re=~A)" reynolds))
        (funcall (intern "XLABEL" vgplot) "u-velocity")
        (funcall (intern "YLABEL" vgplot) "y")
        (funcall format-plot "~A" "set grid")
        (funcall format-plot "~A" "set key top left")
        (funcall format-plot "~A" "set xrange [-0.5:1.1]")
        (funcall format-plot "~A" "set yrange [0:1]")

        ;; Plot reference data (points) and computed data (line)
        (funcall (intern "PLOT" vgplot)
                 ref-u ref-y
                 (format nil "title '~A' with points pt 7 ps 1.5 lc rgb '#E41A1C'" source-name)
                 comp-u comp-y
                 "title 'StreamVorti DCPSE' with lines lw 2 lc rgb '#377EB8'")

        ;; Save plot to file using vgplot:print-plot
        (funcall (intern "PRINT-PLOT" vgplot) png-path)

        (funcall (intern "CLOSE-ALL-PLOTS" vgplot)))

      (format t "Plot saved to: ~A~%" (namestring png-path))
      (namestring png-path))))

;;; ============================================================
;;; ASCII plotting (fallback)
;;; ============================================================

(defun ascii-plot (computed-data reference-data &key (width 60) (height 20))
  "Simple ASCII plot for terminal output.
   Plots u-velocity (x-axis) vs y (y-axis)."
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
    (format t "~%u-velocity vs y (o=reference, *=computed)~%")
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
                         (reynolds 100)
                         (output-dir "plots/")
                         (generate-plots t))
  "Run full validation: load data, compute errors, generate plots.

   Keyword arguments:
     :results-file    Path to StreamVorti centerline data file
     :reynolds        Reynolds number (100, 400, 1000 for Ghia; 1000-21000 for Erturk)
     :output-dir      Directory for output plots (default: plots/)
     :generate-plots  If T, generate PNG plots using vgplot (default: T)"

  ;; Get reference data (may have multiple sources for same Re)
  (multiple-value-bind (reference-list source-names) (load-reference-data reynolds)
    (format t "~%~A~%" (make-string 70 :initial-element #\=))
    (format t "StreamVorti Validation for Re=~A~%" reynolds)
    (format t "~A~%~%" (make-string 70 :initial-element #\=))

    (format t "Reynolds number: ~A~%" reynolds)
    (format t "Reference sources: ~{~A~^, ~}~%" source-names)
    (format t "Results file: ~A~%" results-file)

    ;; Load computed results
    (let ((computed (load-centerline-data results-file)))
      (unless computed
        (format t "~%ERROR: Could not load results file~%")
        (return-from run-validation nil))

      (format t "~%Loaded ~A computed data points~%" (length computed))

      ;; Compute and display errors for each reference source
      (let ((l2-threshold 0.05)
            (linf-threshold 0.10)
            (all-pass t))
        (loop for reference in reference-list
              for source in source-names
              do (format t "~%--- ~A ---~%" source)
                 (format t "Reference data: ~A points~%" (length reference))
                 (multiple-value-bind (l2 linf) (compute-errors computed reference)
                   (format t "  L2 error:   ~,6f~%" l2)
                   (format t "  Linf error: ~,6f~%" linf)
                   (let ((l2-pass (< l2 l2-threshold))
                         (linf-pass (< linf linf-threshold)))
                     (format t "  L2 < ~,0f%: ~A~%" (* 100 l2-threshold) (if l2-pass "PASS" "FAIL"))
                     (format t "  Linf < ~,0f%: ~A~%" (* 100 linf-threshold) (if linf-pass "PASS" "FAIL"))
                     (unless (and l2-pass linf-pass)
                       (setf all-pass nil)))))

        (format t "~%Overall: ~A~%" (if all-pass "PASS" "FAIL")))

      ;; Generate PNG plot if vgplot available
      (when generate-plots
        (generate-plot computed (first reference-list)
                       reynolds (first source-names)
                       :output-dir output-dir))

      ;; Always generate ASCII plot for terminal/logs
      (format t "~%")
      (ascii-plot computed (first reference-list))

      (format t "~%~A~%" (make-string 70 :initial-element #\=))

      ;; Return the computed data for further analysis
      computed)))

;;; Print usage when loaded
(format t "~%StreamVorti Validation Script Loaded~%")
(format t "Usage: (streamvorti.validation:run-validation :results-file \"path/to/centerline.dat\")~%~%")
