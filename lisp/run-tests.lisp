#!/usr/bin/env -S sbcl --script
;;;; run-tests.lisp - Run StreamVorti Lisp test suite
;;;;
;;;; Usage (all equivalent):
;;;;   ./lisp/run-tests.lisp                  # direct (uses shebang)
;;;;   sbcl --script lisp/run-tests.lisp      # explicit SBCL
;;;;   ecl -shell lisp/run-tests.lisp         # explicit ECL
;;;;   (load "lisp/run-tests.lisp")           # from REPL

;;; Determine project directory (where this script lives)
(defvar *script-dir*
  (make-pathname :directory (pathname-directory *load-truename*)))

(defvar *project-dir*
  (make-pathname :directory (butlast (pathname-directory *script-dir*))))

;;; Change to project directory for relative paths to work
(let ((project-path (namestring *project-dir*)))
  #+sbcl (setf *default-pathname-defaults* (pathname project-path))
  #+ecl (ext:chdir project-path))

;;; Load the validation module
(load (merge-pathnames "lisp/compare-ghia.lisp" *project-dir*))

;;; Load the test suite
(load (merge-pathnames "lisp/test-compare-ghia.lisp" *project-dir*))

;;; Exit with appropriate code
(let ((fail-count streamvorti.validation.tests::*fail-count*))
  #+sbcl (sb-ext:exit :code (if (zerop fail-count) 0 1))
  #+ecl (ext:quit (if (zerop fail-count) 0 1))
  #-(or sbcl ecl) (quit))
