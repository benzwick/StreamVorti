;;;; sdl.lisp - SDL v2 entry point
;;;;
;;;; StreamVorti - Software for solving PDEs using explicit methods.
;;;; Copyright (C) 2026 Benjamin F. Zwick
;;;;
;;;; This program is free software: you can redistribute it and/or modify
;;;; it under the terms of the GNU General Public License as published by
;;;; the Free Software Foundation, either version 3 of the License, or
;;;; (at your option) any later version.
;;;;
;;;; Load this file to get the complete SDL v2 API:
;;;;   sbcl --load lisp/sdl.lisp
;;;;   ecl -load lisp/sdl.lisp

(load "lisp/packages.lisp")
(load "lisp/geometry.lisp")
(load "lisp/boundaries.lisp")
(load "lisp/sdl-macros.lisp")
(load "lisp/simulation.lisp")
