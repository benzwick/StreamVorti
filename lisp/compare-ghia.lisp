;;;; compare-ghia.lisp - Validation against benchmark reference data
;;;;
;;;; Copyright (c) 2026 Benjamin F. Zwick
;;;;
;;;; Compares StreamVorti lid-driven cavity results against:
;;;;   - Ghia, Ghia & Shin (1982) for Re <= 1000
;;;;   - Erturk & Corke (2005) for Re > 1000
;;;;
;;;; Usage:
;;;;   sbcl --load ~/quicklisp/setup.lisp \
;;;;        --load lisp/compare-ghia.lisp \
;;;;        --eval '(streamvorti.validation:run-validation :results-file "path/to/data.dat")'

(defpackage :streamvorti.validation
  (:use :cl)
  (:export #:run-validation
           #:load-centerline-data
           #:load-reference-data
           #:compute-errors
           #:generate-plot
           #:ascii-plot))

(in-package :streamvorti.validation)

;;; ============================================================
;;; Configuration
;;; ============================================================

(defparameter *data-dir* "data/")

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

(defun file-ready-p (pathname)
  "Return T if PATHNAME exists and has non-zero size."
  (and (probe-file pathname)
       (ignore-errors
         (plusp (with-open-file (s pathname :element-type '(unsigned-byte 8))
                  (file-length s))))))

(defun await-file (pathname &key (timeout 5.0) (poll-interval 0.1))
  "Wait for PATHNAME to exist with non-zero size.
Returns T if file is ready, NIL on timeout."
  (loop with deadline = (+ (get-internal-real-time)
                           (* timeout internal-time-units-per-second))
        when (file-ready-p pathname) return t
        when (> (get-internal-real-time) deadline) return nil
        do (sleep poll-interval)))

;;; ============================================================
;;; Reference Data
;;; ============================================================

(defun parse-reynolds-header (line)
  "Extract Reynolds numbers from header like '# y Re=100 Re=400'."
  (loop for token in (split-whitespace line)
        when (and (> (length token) 3)
                  (string-equal "Re=" (subseq token 0 3)))
          collect (parse-integer (subseq token 3) :junk-allowed t)))

(defun load-reference-file (filename)
  "Load reference data from FILENAME.
Returns (VALUES reynolds-list data-hash)."
  (let ((path (merge-pathnames filename *data-dir*))
        (reynolds nil)
        (data (make-hash-table)))
    (unless (probe-file path)
      (warn "Reference file not found: ~A" path)
      (return-from load-reference-file (values nil nil)))
    (with-open-file (in path)
      ;; Parse header
      (loop for line = (read-line in nil)
            while (and line (plusp (length line)) (char= #\# (char line 0)))
            do (setf reynolds (or (parse-reynolds-header line) reynolds)))
      ;; Initialize
      (dolist (re reynolds)
        (setf (gethash re data) nil))
      ;; Parse data
      (file-position in 0)
      (loop for line = (read-line in nil)
            while line
            for tokens = (split-whitespace line)
            for y = (parse-number (first tokens))
            when (and y (plusp (length line)) (char/= #\# (char line 0)))
              do (loop for re in reynolds
                       for u in (mapcar #'parse-number (rest tokens))
                       when u do (push (cons y u) (gethash re data)))))
    ;; Restore order
    (maphash (lambda (k v) (setf (gethash k data) (nreverse v))) data)
    (values reynolds data)))

(defvar *ghia-reynolds* nil)
(defvar *ghia-data* nil)
(defvar *erturk-reynolds* nil)
(defvar *erturk-data* nil)

(defun load-reference-data (reynolds)
  "Load reference data for REYNOLDS number.
Returns (VALUES data source-name) or (VALUES NIL NIL)."
  ;; Lazy load
  (unless *ghia-data*
    (multiple-value-setq (*ghia-reynolds* *ghia-data*)
      (load-reference-file "ghia_1982_u.txt")))
  (unless *erturk-data*
    (multiple-value-setq (*erturk-reynolds* *erturk-data*)
      (load-reference-file "erturk_2005_u.txt")))
  ;; Return appropriate source
  (cond
    ((gethash reynolds *ghia-data*)
     (values (gethash reynolds *ghia-data*) "Ghia et al. (1982)"))
    ((gethash reynolds *erturk-data*)
     (values (gethash reynolds *erturk-data*) "Erturk & Corke (2005)"))
    (t
     (warn "No reference data for Re=~A" reynolds)
     (values nil nil))))

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

(defun interpolate-at (y data)
  "Linearly interpolate DATA at Y."
  (let ((sorted (sort (copy-list data) #'< :key #'car)))
    (cond
      ((<= y (caar sorted)) (cdar sorted))
      ((>= y (caar (last sorted))) (cdar (last sorted)))
      (t (loop for (p1 p2) on sorted
               when (and p2 (<= (car p1) y (car p2)))
                 return (let ((t-val (/ (- y (car p1)) (- (car p2) (car p1)))))
                          (+ (cdr p1) (* t-val (- (cdr p2) (cdr p1))))))))))

(defun compute-errors (computed reference)
  "Compute L2 and L-infinity errors. Returns (VALUES l2 linf)."
  (let* ((errors (loop for (y . u-ref) in reference
                       collect (- (interpolate-at y computed) u-ref)))
         (l2 (sqrt (/ (loop for e in errors sum (* e e)) (length errors))))
         (linf (loop for e in errors maximize (abs e))))
    (values l2 linf)))

;;; ============================================================
;;; Plotting with vgplot
;;; ============================================================

(defvar *vgplot-loaded* nil)

(defun ensure-vgplot ()
  "Load vgplot via quicklisp if not already loaded. Returns T on success."
  (when *vgplot-loaded*
    (return-from ensure-vgplot t))
  (handler-case
      (progn
        ;; Load quicklisp if needed
        (unless (find-package :ql)
          (let ((init (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname))))
            (when (probe-file init)
              (load init))))
        ;; Load vgplot
        (when (find-package :ql)
          (funcall (intern "QUICKLOAD" :ql) :vgplot :silent t)
          (setf *vgplot-loaded* t)))
    (error (e)
      (format *error-output* "vgplot not available: ~A~%" e)
      nil)))

(defun vg (name &rest args)
  "Call vgplot function NAME with ARGS."
  (apply (intern (string-upcase (string name)) :vgplot) args))

(defun list-to-vector (list)
  "Convert LIST to a simple vector."
  (coerce list 'vector))

(defun generate-plot (computed reference reynolds source
                      &key (output-dir "plots/") (filename "validation_plot.png"))
  "Generate PNG plot using vgplot. Returns pathname on success, NIL on failure."
  (unless (ensure-vgplot)
    (return-from generate-plot nil))
  (let ((png-path (merge-pathnames filename output-dir)))
    (ensure-directories-exist png-path)
    ;; Extract data as vectors (vgplot expects vectors)
    (let ((ref-u (list-to-vector (mapcar #'cdr reference)))
          (ref-y (list-to-vector (mapcar #'car reference)))
          (comp-u (list-to-vector (mapcar #'cdr computed)))
          (comp-y (list-to-vector (mapcar #'car computed))))
      ;; Plot: reference points and computed line
      (vg 'plot
          ref-u ref-y (format nil "og;~A;" source)
          comp-u comp-y "b;StreamVorti DCPSE;")
      ;; Labels and formatting
      (vg 'title (format nil "Lid-Driven Cavity: u-velocity at x=0.5 (Re=~A)" reynolds))
      (vg 'xlabel "u-velocity")
      (vg 'ylabel "y")
      (vg 'axis '(-0.5 1.1 0 1))
      (vg 'legend :northeast)
      ;; Save to file
      (vg 'print-plot png-path)
      ;; Give gnuplot time to process and write the file
      (sleep 1)
      (vg 'close-all-plots))
    ;; Wait for file to appear (gnuplot may still be writing)
    (if (await-file png-path :timeout 10.0)
        (progn
          (format t "Plot saved: ~A~%" (namestring png-path))
          png-path)
        (progn
          (format *error-output* "Failed to generate plot~%")
          nil))))

;;; ============================================================
;;; ASCII Plot
;;; ============================================================

(defun ascii-plot (computed reference &key (width 60) (height 20))
  "Print ASCII plot to *standard-output*."
  (let* ((all-u (append (mapcar #'cdr computed) (mapcar #'cdr reference)))
         (u-min (reduce #'min all-u))
         (u-max (reduce #'max all-u))
         (u-range (max 1e-6 (- u-max u-min)))
         (canvas (make-array (list height width) :initial-element #\Space)))
    ;; Reference as 'o'
    (loop for (y . u) in reference
          for row = (min (1- height) (max 0 (floor (* (- 1.0 y) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- u u-min) u-range) (1- width)))))
          do (setf (aref canvas row col) #\o))
    ;; Computed as '*'
    (loop for (y . u) in computed
          for row = (min (1- height) (max 0 (floor (* (- 1.0 y) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- u u-min) u-range) (1- width)))))
          do (setf (aref canvas row col) #\*))
    ;; Print
    (format t "~%u-velocity vs y (o=reference, *=computed)~%")
    (format t "u: ~,2F~VT~,2F~%" u-min width u-max)
    (dotimes (row height)
      (format t "|")
      (dotimes (col width)
        (write-char (aref canvas row col)))
      (format t "|~%"))
    (format t "y: 1.0~VT0.0~%" (- width 3))))

;;; ============================================================
;;; Main Entry Point
;;; ============================================================

(defun run-validation (&key
                         (results-file "output_dat/mfem_square10x10_u_centerline_x0.5.dat")
                         (reynolds 100)
                         (output-dir "plots/")
                         (plot t))
  "Validate StreamVorti results against reference data.
Returns T if validation passes, NIL otherwise."
  (format t "~%~70,,,'-<~>~%")
  (format t "StreamVorti Validation for Re=~A~%" reynolds)
  (format t "~70,,,'-<~>~%~%")
  (format t "Reynolds number: ~A~%" reynolds)
  (format t "Results file: ~A~%" results-file)
  ;; Load data
  (let ((computed (load-centerline-data results-file)))
    (format t "~%Loaded ~A computed data points~%" (length computed))
    (multiple-value-bind (reference source) (load-reference-data reynolds)
      (unless reference
        (format t "~%ERROR: No reference data for Re=~A~%" reynolds)
        (return-from run-validation nil))
      (format t "~%--- ~A ---~%" source)
      (format t "Reference data: ~A points~%" (length reference))
      ;; Compute errors
      (multiple-value-bind (l2 linf) (compute-errors computed reference)
        (let ((l2-pass (< l2 0.05))
              (linf-pass (< linf 0.10)))
          (format t "  L2 error:   ~,6F~%" l2)
          (format t "  Linf error: ~,6F~%" linf)
          (format t "  L2 < 5%: ~A~%" (if l2-pass "PASS" "FAIL"))
          (format t "  Linf < 10%: ~A~%" (if linf-pass "PASS" "FAIL"))
          (format t "~%Overall: ~A~%" (if (and l2-pass linf-pass) "PASS" "FAIL"))
          ;; Plots
          (when plot
            (generate-plot computed reference reynolds source :output-dir output-dir))
          (ascii-plot computed reference)
          (format t "~%~70,,,'-<~>~%")
          (and l2-pass linf-pass))))))

;;; ============================================================

(format t "~%StreamVorti Validation Script Loaded~%")
(format t "Usage: (streamvorti.validation:run-validation :results-file \"path/to/centerline.dat\")~%~%")
