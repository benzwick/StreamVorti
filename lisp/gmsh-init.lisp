;;;; gmsh-init.lisp - Bootstrap gmsh-cl in embedded ECL
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.
;;;;
;;;; Loaded conditionally when STREAMVORTI_WITH_GMSH is enabled.
;;;; Expects *gmsh-cl-dir* to be set before loading this file.

(in-package :cl-user)

(let ((gmsh-dir (truename *gmsh-cl-dir*)))

  ;; Register gmsh-cl source with ASDF
  (push gmsh-dir asdf:*central-registry*)

  ;; Register ocicl dependencies (cffi, alexandria, etc.)
  ;; ocicl install places each dependency in a subdirectory of ocicl/
  (let ((ocicl-dir (merge-pathnames "ocicl/" gmsh-dir)))
    (when (probe-file ocicl-dir)
      (dolist (entry (directory (merge-pathnames "*/" ocicl-dir)))
        (push entry asdf:*central-registry*))))

  ;; Load CFFI first so we can set library search paths
  (asdf:load-system :cffi)

  ;; Tell CFFI where to find libgmsh.so (built from gmsh-cl's submodule)
  (let ((gmsh-lib-dir (merge-pathnames "_reference/gmsh/build/" gmsh-dir)))
    (when (probe-file gmsh-lib-dir)
      (pushnew gmsh-lib-dir
               (symbol-value (intern "*FOREIGN-LIBRARY-DIRECTORIES*" "CFFI"))
               :test #'equal)))

  ;; Load gmsh-cl and its dependencies.
  ;; gmsh-cl's defpackage forms will update the stub packages created
  ;; by packages.lisp, adding exports and full functionality.
  (asdf:load-system :gmsh-cl))
