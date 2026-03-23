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

;;; Load the validation module and tests
(load (merge-pathnames "lisp/compare-ghia.lisp" *project-dir*))
(load (merge-pathnames "lisp/test-compare-ghia.lisp" *project-dir*))

;;; Load SDL packages (needed for test-sdl.lisp)
(load (merge-pathnames "lisp/packages.lisp" *project-dir*))
(load (merge-pathnames "lisp/geometry.lisp" *project-dir*))
(load (merge-pathnames "lisp/mesh.lisp" *project-dir*))
(load (merge-pathnames "lisp/boundaries.lisp" *project-dir*))
(load (merge-pathnames "lisp/sdl-macros.lisp" *project-dir*))
(load (merge-pathnames "lisp/simulation.lisp" *project-dir*))
(load (merge-pathnames "lisp/cpp-bridge.lisp" *project-dir*))

;;; Load SDL test suite
(load (merge-pathnames "lisp/test-sdl.lisp" *project-dir*))

;;; Exit with appropriate code
(let ((total-fails (+ streamvorti.validation.tests::*fail-count*
                      streamvorti.sdl.tests::*fail-count*)))
  #+sbcl (sb-ext:exit :code (if (zerop total-fails) 0 1))
  #+ecl (ext:quit (if (zerop total-fails) 0 1))
  #-(or sbcl ecl) (quit))
