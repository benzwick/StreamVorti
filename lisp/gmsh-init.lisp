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

;; Delete stub packages created by packages.lisp so gmsh-cl can
;; define the real ones with full exports and nicknames.
(dolist (name '(:gmsh/fltk :gmsh/option :gmsh/mesh :gmsh/geo
               :gmsh/occ :gmsh))
  (let ((pkg (find-package name)))
    (when (and pkg (null (apropos-list "" pkg)))
      (delete-package pkg))))

;; Register gmsh-cl source with ASDF
(push (truename *gmsh-cl-dir*) asdf:*central-registry*)

;; Load gmsh-cl and its dependencies
(asdf:load-system :gmsh-cl)
