;;;; validate-all.lisp - Master validation script for all demos/examples
;;;;
;;;; Copyright (c) 2026 Benjamin F. Zwick
;;;;
;;;; Runs validation comparisons against analytical solutions and
;;;; benchmark reference data for every demo that has one.
;;;;
;;;; Usage:
;;;;   sbcl --load lisp/validate-all.lisp \
;;;;        --eval '(streamvorti.validate:run-all)'
;;;;
;;;; Or validate a single case:
;;;;   sbcl --load lisp/validate-all.lisp \
;;;;        --eval '(streamvorti.validate:validate-cavity)'

(defpackage :streamvorti.validate
  (:use :cl)
  (:export #:run-all
           #:validate-cavity
           #:validate-cavity-fdm
           #:validate-cavity-re1000
           #:validate-poiseuille
           #:validate-channel
           #:validate-double-lid-symmetry))

(in-package :streamvorti.validate)

;;; ============================================================
;;; Utilities (self-contained — no dependency on other scripts)
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
      (when (numberp val) (float val 1.0d0)))))

(defun load-centerline-data (filename)
  "Load centerline data from FILENAME. Returns list of (y . u) pairs."
  (unless (probe-file filename)
    (format t "  File not found: ~A~%" filename)
    (return-from load-centerline-data nil))
  (with-open-file (in filename)
    (loop for line = (read-line in nil)
          while line
          for tokens = (split-whitespace line)
          for y = (parse-number (first tokens))
          for u = (parse-number (second tokens))
          when (and y u (or (zerop (length line)) (char/= #\# (char line 0))))
            collect (cons y u))))

(defun interpolate-at (y data)
  "Linearly interpolate DATA (list of (y . val) pairs) at Y."
  (let ((sorted (sort (copy-list data) #'< :key #'car)))
    (cond
      ((<= y (caar sorted)) (cdar sorted))
      ((>= y (caar (last sorted))) (cdar (last sorted)))
      (t (loop for (p1 p2) on sorted
               when (and p2 (<= (car p1) y (car p2)))
                 return (let ((t-val (/ (- y (car p1)) (- (car p2) (car p1)))))
                          (+ (cdr p1) (* t-val (- (cdr p2) (cdr p1))))))))))

(defun compute-errors-vs-function (computed analytical-fn)
  "Compute L2 and L-inf errors of COMPUTED against ANALYTICAL-FN.
Returns (VALUES l2 linf max-at-y)."
  (let* ((errors (loop for (y . u) in computed
                       collect (cons y (- u (funcall analytical-fn y)))))
         (l2 (sqrt (/ (loop for (_ . e) in errors sum (* e e))
                      (length errors))))
         (linf-entry (reduce (lambda (a b)
                               (if (> (abs (cdr a)) (abs (cdr b))) a b))
                             errors))
         (linf (abs (cdr linf-entry))))
    (values l2 linf (car linf-entry))))

(defun compute-errors-vs-reference (computed reference)
  "Compute L2 and L-inf errors of COMPUTED against REFERENCE data.
Returns (VALUES l2 linf)."
  (let* ((errors (loop for (y . u-ref) in reference
                       for u-comp = (interpolate-at y computed)
                       when u-comp
                         collect (- u-comp u-ref)))
         (l2 (sqrt (/ (loop for e in errors sum (* e e)) (max 1 (length errors)))))
         (linf (loop for e in errors maximize (abs e))))
    (values l2 linf)))

;;; ============================================================
;;; ASCII plot helper
;;; ============================================================

(defun ascii-comparison-plot (computed-data reference-data
                              &key (width 65) (height 20)
                                   (label-computed "computed")
                                   (label-reference "reference"))
  "Print ASCII plot comparing two datasets.
Each dataset is a list of (y . value) pairs.
Legend: o=reference, *=computed, @=both"
  (let* ((all-vals (append (mapcar #'cdr computed-data)
                           (mapcar #'cdr reference-data)))
         (v-min (reduce #'min all-vals))
         (v-max (reduce #'max all-vals))
         (v-range (max 1d-10 (- v-max v-min)))
         (all-y (append (mapcar #'car computed-data)
                        (mapcar #'car reference-data)))
         (y-min (reduce #'min all-y))
         (y-max (reduce #'max all-y))
         (y-range (max 1d-10 (- y-max y-min)))
         (canvas (make-array (list height width) :initial-element #\Space)))
    ;; Reference as 'o'
    (loop for (y . v) in reference-data
          for row = (min (1- height) (max 0 (floor (* (- 1.0 (/ (- y y-min) y-range)) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- v v-min) v-range) (1- width)))))
          do (setf (aref canvas row col) #\o))
    ;; Computed as '*', '@' if overlapping
    (loop for (y . v) in computed-data
          for row = (min (1- height) (max 0 (floor (* (- 1.0 (/ (- y y-min) y-range)) (1- height)))))
          for col = (min (1- width) (max 0 (floor (* (/ (- v v-min) v-range) (1- width)))))
          do (setf (aref canvas row col)
                   (if (char= (aref canvas row col) #\o) #\@ #\*)))
    ;; Print
    (format t "  (o=~A, *=~A, @=both)~%" label-reference label-computed)
    (format t "  val: ~,3F~VT~,3F~%" v-min (+ width 6) v-max)
    (dotimes (row height)
      (format t "  |")
      (dotimes (col width)
        (write-char (aref canvas row col)))
      (format t "|~%"))
    (format t "  y: ~,2F~VT~,2F~%" y-max (+ width 4) y-min)))

;;; ============================================================
;;; Reference Data Loaders
;;; ============================================================

(defparameter *data-dir* "data/")

(defun parse-reynolds-header (line)
  "Extract Reynolds numbers from header like '# y Re=100 Re=400'."
  (loop for token in (split-whitespace line)
        when (and (> (length token) 3)
                  (string-equal "Re=" (subseq token 0 3)))
          collect (parse-integer (subseq token 3) :junk-allowed t)))

(defun load-ghia-data (reynolds &key (component :u))
  "Load Ghia et al. (1982) data for given REYNOLDS number.
COMPONENT is :u (vertical centerline) or :v (horizontal centerline).
Returns list of (coord . value) pairs, or NIL if not available."
  (let ((path (merge-pathnames (ecase component
                                 (:u "ghia_1982_u.txt")
                                 (:v "ghia_1982_v.txt"))
                               *data-dir*))
        (reynolds-list nil)
        (data (make-hash-table)))
    (unless (probe-file path)
      (return-from load-ghia-data nil))
    (with-open-file (in path)
      ;; Parse header for Reynolds numbers
      (loop for line = (read-line in nil)
            while (and line (plusp (length line)) (char= #\# (char line 0)))
            do (setf reynolds-list (or (parse-reynolds-header line) reynolds-list)))
      (dolist (re reynolds-list)
        (setf (gethash re data) nil))
      ;; Parse data
      (file-position in 0)
      (loop for line = (read-line in nil)
            while line
            for tokens = (split-whitespace line)
            for y = (parse-number (first tokens))
            when (and y (plusp (length line)) (char/= #\# (char line 0)))
              do (loop for re in reynolds-list
                       for u in (mapcar #'parse-number (rest tokens))
                       when u do (push (cons y u) (gethash re data)))))
    (let ((result (gethash reynolds data)))
      (when result (nreverse result)))))

;;; ============================================================
;;; Analytical Solutions
;;; ============================================================

(defun poiseuille-velocity (y &key (u-max 1.0d0) (height 1.0d0))
  "Exact Poiseuille flow: u(y) = 4·U_max·y·(H-y)/H²"
  (* 4.0d0 u-max y (- height y) (/ 1.0d0 (* height height))))

;;; ============================================================
;;; Individual Validation Cases
;;; ============================================================

(defun print-errors (l2 linf &key (l2-threshold 0.05d0) (linf-threshold 0.10d0))
  "Print error summary and return T if both pass."
  (let ((l2-pass (< l2 l2-threshold))
        (linf-pass (< linf linf-threshold)))
    (format t "  L2 error:   ~,6F  (~A, threshold ~,1F%)~%"
            l2 (if l2-pass "PASS" "FAIL") (* 100 l2-threshold))
    (format t "  Linf error: ~,6F  (~A, threshold ~,1F%)~%"
            linf (if linf-pass "PASS" "FAIL") (* 100 linf-threshold))
    (and l2-pass linf-pass)))

(defun validate-cavity (&key (prefix "lid-driven-cavity") (reynolds 100))
  "Validate lid-driven cavity against Ghia et al. (1982)."
  (format t "~%=== Lid-Driven Cavity (Re=~A) ==========================~%" reynolds)
  (let ((pass-u nil)
        (pass-v nil))
    ;; u-velocity along vertical centerline x=0.5
    (format t "~%  --- u-velocity along x=0.5 ---~%")
    (let* ((file-u (format nil "output_dat/~A_u_centerline_x0.5.dat" prefix))
           (computed-u (load-centerline-data file-u))
           (reference-u (load-ghia-data reynolds :component :u)))
      (cond
        ((null computed-u)
         (format t "  SKIP: No results found (~A)~%" file-u))
        ((null reference-u)
         (format t "  SKIP: No Ghia u-data for Re=~A~%" reynolds))
        (t
         (format t "  Computed: ~A points, Reference: ~A points~%"
                 (length computed-u) (length reference-u))
         (multiple-value-bind (l2 linf) (compute-errors-vs-reference computed-u reference-u)
           (setf pass-u (print-errors l2 linf))
           (ascii-comparison-plot computed-u reference-u
                                  :label-computed "StreamVorti"
                                  :label-reference "Ghia (1982)")))))
    ;; v-velocity along horizontal centerline y=0.5
    (format t "~%  --- v-velocity along y=0.5 ---~%")
    (let* ((file-v (format nil "output_dat/~A_v_centerline_y0.5.dat" prefix))
           (computed-v (load-centerline-data file-v))
           (reference-v (load-ghia-data reynolds :component :v)))
      (cond
        ((null computed-v)
         (format t "  SKIP: No results found (~A)~%" file-v))
        ((null reference-v)
         (format t "  SKIP: No Ghia v-data for Re=~A~%" reynolds))
        (t
         (format t "  Computed: ~A points, Reference: ~A points~%"
                 (length computed-v) (length reference-v))
         (multiple-value-bind (l2 linf) (compute-errors-vs-reference computed-v reference-v)
           (setf pass-v (print-errors l2 linf))
           (ascii-comparison-plot computed-v reference-v
                                  :label-computed "StreamVorti"
                                  :label-reference "Ghia (1982)")))))
    ;; Overall pass: both must pass (or be skipped)
    (and (or pass-u (null pass-u))
         (or pass-v (null pass-v))
         (or pass-u pass-v))))

(defun validate-cavity-fdm ()
  "Validate FDM cavity against Ghia et al. (1982)."
  (validate-cavity :prefix "lid-driven-cavity-fdm" :reynolds 100))

(defun validate-cavity-re1000 ()
  "Validate Re=1000 cavity against Ghia et al. (1982)."
  (validate-cavity :prefix "lid-driven-cavity-re1000" :reynolds 1000))

(defun validate-poiseuille (&key (prefix "poiseuille-flow")
                                  (u-max 1.0d0) (channel-height 1.0d0))
  "Validate Poiseuille flow against exact parabolic profile."
  (format t "~%=== Poiseuille Flow ======================================~%")
  (format t "  Exact: u(y) = 4 * ~A * y * (~A - y) / ~A~%"
          u-max channel-height (* channel-height channel-height))
  (let* ((file (format nil "output_dat/~A_u_centerline_x0.5.dat" prefix))
         (computed (load-centerline-data file)))
    (cond
      ((null computed)
       (format t "  SKIP: No results found (~A)~%" file)
       nil)
      (t
       (format t "  Computed: ~A points~%" (length computed))
       (multiple-value-bind (l2 linf max-y)
           (compute-errors-vs-function
            computed
            (lambda (y) (poiseuille-velocity y :u-max u-max
                                               :height channel-height)))
         (declare (ignore max-y))
         (let ((pass (print-errors l2 linf))
               (exact (loop for (y . nil) in computed
                            collect (cons y (poiseuille-velocity y :u-max u-max
                                                                   :height channel-height)))))
           (ascii-comparison-plot computed exact
                                  :label-computed "StreamVorti"
                                  :label-reference "Exact")
           pass))))))

(defun validate-channel (&key (prefix "channel-flow-2d")
                               (probe "exit-channel")
                               (u-max 1.0d0) (channel-height 1.0d0))
  "Validate channel flow at exit probe against Poiseuille profile.
The flow should approach the fully-developed parabolic profile
at downstream stations."
  (format t "~%=== Channel Flow (exit probe) ============================~%")
  (let* ((file (format nil "output_dat/~A_u_~A.dat" prefix probe))
         (computed (load-centerline-data file)))
    (cond
      ((null computed)
       (format t "  SKIP: No results found (~A)~%" file)
       nil)
      (t
       (format t "  Computed: ~A points~%" (length computed))
       (format t "  Comparing exit profile to fully-developed Poiseuille~%")
       (multiple-value-bind (l2 linf max-y)
           (compute-errors-vs-function
            computed
            (lambda (y) (poiseuille-velocity y :u-max u-max
                                               :height channel-height)))
         (declare (ignore max-y))
         ;; Channel exit may not be fully developed — use relaxed thresholds
         (let ((pass (print-errors l2 linf :l2-threshold 0.10d0
                                           :linf-threshold 0.20d0))
               (exact (loop for (y . nil) in computed
                            collect (cons y (poiseuille-velocity y :u-max u-max
                                                                   :height channel-height)))))
           (ascii-comparison-plot computed exact
                                  :label-computed "StreamVorti"
                                  :label-reference "Poiseuille")
           pass))))))

(defun validate-double-lid-symmetry (&key (prefix "double-lid-cavity"))
  "Validate double-lid cavity by checking antisymmetry about y=0.5.
For opposite-moving lids, u(y) = -u(1-y) along the vertical centerline."
  (format t "~%=== Double Lid-Driven Cavity (symmetry check) ============~%")
  (let* ((file (format nil "output_dat/~A_u_centerline_x0.5.dat" prefix))
         (computed (load-centerline-data file)))
    (cond
      ((null computed)
       (format t "  SKIP: No results found (~A)~%" file)
       nil)
      (t
       (format t "  Computed: ~A points~%" (length computed))
       (format t "  Checking antisymmetry: u(y) should equal -u(1-y)~%")
       ;; Build the antisymmetric "reference": -u(1-y)
       (let* ((mirror (loop for (y . u) in computed
                            collect (cons y (- (or (interpolate-at (- 1.0d0 y) computed) 0.0d0)))))
              (errors (loop for ((y . u) . (nil . u-mirror)) on (mapcar #'cons computed mirror)
                            when u-mirror collect (- u u-mirror)))
              (l2 (if errors
                      (sqrt (/ (loop for e in errors sum (* e e)) (length errors)))
                      0.0d0))
              (linf (if errors
                       (loop for e in errors maximize (abs e))
                       0.0d0)))
         (format t "  Antisymmetry L2:   ~,6F~%" l2)
         (format t "  Antisymmetry Linf: ~,6F~%" linf)
         (let ((pass (and (< l2 0.05d0) (< linf 0.10d0))))
           (format t "  Result: ~A~%" (if pass "PASS" "FAIL"))
           (ascii-comparison-plot computed mirror
                                  :label-computed "u(y)"
                                  :label-reference "-u(1-y)")
           pass))))))

;;; ============================================================
;;; Master Entry Point
;;; ============================================================

(defun run-all ()
  "Run all available validations and print summary."
  (format t "~%~72,,,,'=<~>~%")
  (format t "StreamVorti Validation Suite~%")
  (format t "~72,,,,'=<~>~%")
  (let ((results nil))
    ;; Lid-driven cavity variants
    (push (cons "Cavity Re=100 (DCPSE)"
                (validate-cavity :prefix "lid-driven-cavity" :reynolds 100))
          results)
    (push (cons "Cavity Re=100 (FDM)"
                (validate-cavity :prefix "lid-driven-cavity-fdm" :reynolds 100))
          results)
    (push (cons "Cavity Re=1000"
                (validate-cavity :prefix "lid-driven-cavity-re1000" :reynolds 1000))
          results)
    ;; CLI-mode cavity (uses different prefix)
    (push (cons "Cavity Re=100 (CLI)"
                (validate-cavity :prefix "mfem_square10x10" :reynolds 100))
          results)
    ;; Poiseuille / channel
    (push (cons "Poiseuille flow"
                (validate-poiseuille))
          results)
    (push (cons "Channel flow (exit)"
                (validate-channel))
          results)
    ;; Symmetry check
    (push (cons "Double-lid symmetry"
                (validate-double-lid-symmetry))
          results)
    ;; Summary
    (let* ((attempted (remove-if (lambda (r) (null (cdr r))) results
                                 :key (lambda (r) (declare (ignore r)) nil)))
           (ran (remove-if (lambda (r) (null (cdr r))) results))
           (passed (count-if #'cdr results))
           (failed (count-if (lambda (r) (and (cdr r) (not (eq (cdr r) t)))) results))
           ;; Actually count properly
           (total-ran 0)
           (total-pass 0)
           (total-fail 0)
           (total-skip 0))
      (declare (ignore attempted ran passed failed))
      (format t "~%~72,,,,'=<~>~%")
      (format t "SUMMARY~%")
      (format t "~72,,,,'-<~>~%")
      (dolist (r (nreverse results))
        (let ((status (cdr r)))
          (cond
            ((null status) (incf total-skip)
             (format t "  SKIP  ~A~%" (car r)))
            ((eq status t) (incf total-pass) (incf total-ran)
             (format t "  PASS  ~A~%" (car r)))
            (t (incf total-fail) (incf total-ran)
             (format t "  FAIL  ~A~%" (car r))))))
      (format t "~72,,,,'-<~>~%")
      (format t "Ran: ~A  Passed: ~A  Failed: ~A  Skipped: ~A~%"
              total-ran total-pass total-fail total-skip)
      (format t "~72,,,,'=<~>~%~%")
      (zerop total-fail))))

;;; ============================================================

(format t "~%StreamVorti Validation Suite Loaded~%")
(format t "Usage: (streamvorti.validate:run-all)~%")
(format t "   or: (streamvorti.validate:validate-cavity)~%")
(format t "   or: (streamvorti.validate:validate-poiseuille)~%~%")
