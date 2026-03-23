;;;; plot-sparsity.lisp - ASCII sparsity pattern plots via vgplot/gnuplot
;;;;
;;;; Copyright (c) 2026 Benjamin F. Zwick
;;;;
;;;; Reads a MATLAB-format sparse matrix file (as produced by
;;;; MFEM's SparseMatrix::PrintMatlab) and plots the sparsity pattern
;;;; using gnuplot's "dumb" terminal for ASCII art output.
;;;;
;;;; The file format is one triplet per line:  row  col  value
;;;; (1-indexed, as per MATLAB convention)
;;;;
;;;; Usage:
;;;;   sbcl --load ~/quicklisp/setup.lisp \
;;;;        --load lisp/plot-sparsity.lisp \
;;;;        --eval '(streamvorti.sparsity:plot-sparsity-pattern "output_dat/dx.dat")'
;;;;
;;;; Or from a running Lisp session:
;;;;   (streamvorti.sparsity:plot-sparsity-pattern "output_dat/dx.dat")
;;;;   (streamvorti.sparsity:plot-sparsity-pattern "output_dat/dx.dat" :zoom 10)

(defpackage :streamvorti.sparsity
  (:use :cl)
  (:export #:load-sparse-matrix
           #:ascii-sparsity-pattern
           #:plot-sparsity-pattern))

(in-package :streamvorti.sparsity)

;;; ============================================================
;;; Utilities
;;; ============================================================

(defun split-whitespace (string)
  "Split STRING on whitespace, returning list of non-empty tokens."
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
      (when (numberp val) val))))

;;; ============================================================
;;; Sparse Matrix Loading
;;; ============================================================

(defun load-sparse-matrix (filename)
  "Load a MATLAB-format sparse matrix from FILENAME.
Returns (VALUES rows cols nrows ncols nnz) where rows and cols are
vectors of 0-indexed row/column indices of non-zero entries."
  (unless (probe-file filename)
    (error "Sparse matrix file not found: ~A" filename))
  (let ((row-list nil)
        (col-list nil)
        (max-row 0)
        (max-col 0)
        (count 0))
    (with-open-file (in filename)
      (loop for line = (read-line in nil)
            while line
            for tokens = (split-whitespace line)
            when (>= (length tokens) 3)
              do (let ((r (parse-number (first tokens)))
                       (c (parse-number (second tokens))))
                   (when (and r c)
                     (let ((ri (1- (truncate r)))   ; convert to 0-indexed
                           (ci (1- (truncate c))))
                       (push ri row-list)
                       (push ci col-list)
                       (setf max-row (max max-row ri))
                       (setf max-col (max max-col ci))
                       (incf count))))))
    (values (coerce (nreverse row-list) 'vector)
            (coerce (nreverse col-list) 'vector)
            (1+ max-row)
            (1+ max-col)
            count)))

;;; ============================================================
;;; Pure ASCII Sparsity Pattern (no gnuplot dependency)
;;; ============================================================

(defun ascii-sparsity-pattern (rows cols nrows ncols
                               &key (width 72) (height 36)
                                    (title "Sparsity Pattern")
                                    (max-row nil) (max-col nil)
                                    (stream *standard-output*))
  "Print an ASCII sparsity pattern to STREAM.
ROWS and COLS are vectors of 0-indexed non-zero positions.
MAX-ROW/MAX-COL limit the view to a submatrix (for zooming)."
  (let* ((view-rows (or max-row nrows))
         (view-cols (or max-col ncols))
         (nnz (length rows))
         (canvas (make-array (list height width) :initial-element #\Space))
         (plotted 0))
    ;; Map non-zero entries to canvas
    (dotimes (k nnz)
      (let ((r (aref rows k))
            (c (aref cols k)))
        (when (and (< r view-rows) (< c view-cols))
          (let ((cy (min (1- height) (floor (* r (/ height view-rows)))))
                (cx (min (1- width)  (floor (* c (/ width  view-cols))))))
            (setf (aref canvas cy cx) #\.)
            (incf plotted)))))
    ;; Header
    (format stream "~%~A~%" title)
    (format stream "Matrix: ~Dx~D, nnz=~D" nrows ncols nnz)
    (when (or max-row max-col)
      (format stream " (showing ~Dx~D)" view-rows view-cols))
    (format stream "~%")
    ;; Column axis
    (format stream "     0")
    (format stream "~VT~D~%" (+ width 4) (1- view-cols))
    ;; Canvas
    (dotimes (row height)
      (let ((matrix-row (floor (* row (/ view-rows height)))))
        (format stream "~4D " matrix-row))
      (write-char #\| stream)
      (dotimes (col width)
        (write-char (aref canvas row col) stream))
      (write-char #\| stream)
      (terpri stream))
    ;; Bottom axis
    (format stream "~4D " (1- view-rows))
    (format stream "+")
    (dotimes (i width) (write-char #\- stream))
    (format stream "+~%")
    (format stream "     Plotted ~D/~D non-zeros in view~%~%" plotted nnz)))

;;; ============================================================
;;; vgplot-based plotting (gnuplot dumb terminal)
;;; ============================================================

(defvar *vgplot-loaded* nil)

(defun ensure-vgplot ()
  "Load vgplot via quicklisp if not already loaded. Returns T on success."
  (when *vgplot-loaded*
    (return-from ensure-vgplot t))
  (handler-case
      (progn
        (unless (find-package :ql)
          (let ((init (merge-pathnames "quicklisp/setup.lisp"
                                       (user-homedir-pathname))))
            (when (probe-file init)
              (load init))))
        (when (find-package :ql)
          (funcall (intern "QUICKLOAD" :ql) :vgplot :silent t)
          (setf *vgplot-loaded* t)))
    (error (e)
      (format *error-output* "vgplot not available: ~A~%" e)
      nil)))

(defun vg (name &rest args)
  "Call vgplot function NAME with ARGS."
  (apply (intern (string-upcase (string name)) :vgplot) args))

(defun plot-with-vgplot (rows cols nrows ncols
                         &key (title "Sparsity Pattern")
                              (max-row nil) (max-col nil)
                              (width 80) (height 40))
  "Plot sparsity pattern using vgplot with gnuplot's dumb terminal."
  (unless (ensure-vgplot)
    (format t "vgplot not available, falling back to built-in ASCII plot~%")
    (return-from plot-with-vgplot nil))
  (let* ((view-rows (or max-row nrows))
         (view-cols (or max-col ncols))
         (filtered-c nil)
         (filtered-r nil))
    ;; Filter entries within view
    (dotimes (k (length rows))
      (let ((r (aref rows k))
            (c (aref cols k)))
        (when (and (< r view-rows) (< c view-cols))
          (push c filtered-c)
          ;; Flip row for plotting (row 0 at top)
          (push (- view-rows r 1) filtered-r))))
    (let ((cv (coerce (nreverse filtered-c) 'vector))
          (rv (coerce (nreverse filtered-r) 'vector)))
      ;; Set dumb terminal for ASCII output
      (vg 'format-plot nil
          (format nil "set terminal dumb ~D ~D" width height))
      (vg 'format-plot nil "set size ratio -1")
      (vg 'plot cv rv
          (format nil "+k;~A (~Dx~D, nnz=~D);"
                  title view-rows view-cols (length cv)))
      (vg 'xlabel "column")
      (vg 'ylabel "row")
      (vg 'axis (list 0 (1- view-cols)
                      0 (1- view-rows)))
      (vg 'close-all-plots))
    t))

;;; ============================================================
;;; Main Entry Point
;;; ============================================================

(defun plot-sparsity-pattern (filename &key (zoom nil)
                                            (width 72)
                                            (height 36)
                                            (use-vgplot t))
  "Plot the sparsity pattern of a MATLAB-format sparse matrix file.

FILENAME  Path to the matrix file (MFEM PrintMatlab format).
ZOOM      If given, show only the top-left ZOOM x ZOOM submatrix.
WIDTH     ASCII canvas width in characters (default 72).
HEIGHT    ASCII canvas height in lines (default 36).
USE-VGPLOT  If T, try vgplot first; fall back to built-in ASCII on failure.

Examples:
  ;; Full sparsity pattern
  (plot-sparsity-pattern \"output_dat/dx.dat\")

  ;; Zoomed to top-left 10x10 corner
  (plot-sparsity-pattern \"output_dat/dx.dat\" :zoom 10)

  ;; Both full and zoomed
  (plot-sparsity-pattern \"output_dat/dx.dat\" :zoom 15 :width 80 :height 40)"
  (format t "~%Loading sparse matrix from ~A...~%" filename)
  (multiple-value-bind (rows cols nrows ncols nnz)
      (load-sparse-matrix filename)
    (format t "Matrix size: ~D x ~D, non-zeros: ~D~%" nrows ncols nnz)
    (format t "Density: ~,4F%~%" (* 100.0 (/ nnz (* nrows ncols))))
    ;; Full pattern
    (let ((title (format nil "~A - Full View" (file-namestring filename))))
      (if (and use-vgplot (ensure-vgplot))
          (plot-with-vgplot rows cols nrows ncols :title title
                            :width width :height height)
          (ascii-sparsity-pattern rows cols nrows ncols
                                  :title title
                                  :width width :height height)))
    ;; Zoomed view
    (when zoom
      (let ((z (min zoom nrows ncols))
            (title-z (format nil "~A - Top-Left ~Dx~D"
                             (file-namestring filename) zoom zoom)))
        (if (and use-vgplot (ensure-vgplot))
            (plot-with-vgplot rows cols nrows ncols
                              :title title-z
                              :max-row z :max-col z
                              :width width :height height)
            (ascii-sparsity-pattern rows cols nrows ncols
                                    :title title-z
                                    :max-row z :max-col z
                                    :width width :height height))))))

;;; ============================================================

(format t "~%StreamVorti Sparsity Pattern Plotter Loaded~%")
(format t "Usage: (streamvorti.sparsity:plot-sparsity-pattern \"path/to/matrix.dat\")~%")
(format t "       (streamvorti.sparsity:plot-sparsity-pattern \"path/to/matrix.dat\" :zoom 10)~%~%")
