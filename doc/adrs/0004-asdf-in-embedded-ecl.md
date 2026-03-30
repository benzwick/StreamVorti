# ADR-0004: ASDF in Embedded ECL

**Status:** Proposed
**Date:** 2026-03-23
**Authors:** Benjamin F. Zwick

---

## Context

The embedded ECL runtime loads StreamVorti's Lisp files directly via
`(load "file.lisp")`. External Common Lisp libraries like gmsh-cl are
distributed as ASDF systems with dependency management. Loading them requires
ASDF to be available in the runtime.

## Decision

Enable ASDF unconditionally in the embedded ECL by calling `(require 'asdf)`
during `Runtime::init()`, after `cl_boot()` and before loading core SDL files.

ECL ships ASDF as a contributed module — no external installation is required.
The call is unconditional (not gated behind `STREAMVORTI_WITH_GMSH`) because
ASDF is lightweight and may be useful for future library integrations.

## Consequences

- **Slight startup overhead.** Loading ASDF adds a small amount of time to
  ECL initialization. This is negligible compared to the overall solver runtime.
- **Enables loading any ASDF system.** Future CL libraries can be integrated
  without additional infrastructure changes.
- **No impact on existing functionality.** ASDF does not interfere with the
  direct `(load ...)` approach used for core SDL files.
