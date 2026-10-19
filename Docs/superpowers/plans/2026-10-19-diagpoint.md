# DiagPoint implementation plan

**Goal:** Add a runtime fixed-position cell sampler on branch `codex/diagpoint`, based on `development`.

**Approved design:** The conversation specifies `position`, `field_names`, `int`, and `file`; select the finest valid cell and append a time history. Support 2D and 3D, MPI ownership, regridding, and restart. No spatial interpolation.

**Architecture:** Following the user's follow-up, implement `Source/Utility/Diagnostics/DiagPoint.H` and `.cpp` in PelePhysics as a registered `DiagBase` subclass. Include it in PelePhysics's GNU Make package and PeleC's CMake physics library. `volume_fraction_field` is optional for host portability; selecting `vfrac` enables PeleC covered-cell rejection. Record requested coordinates, cell centers, level, step, time, and requested fields. Reject unsupported time-period scheduling and filters. Check the output schema/position on append; skip an identical final step/time and reject older records. Shared tests live in PelePhysics `Testing/Diagnostics/DiagPoint`; the production smoke test remains in PeleC.

**Tech stack:** C++17-compatible AMReX APIs; existing PeleC C++20 build; standalone CMake/CTest tests linked to the checked-out AMReX and Diagnostics base.

- [x] Write a standalone AMReX behavior test before implementation. Verify that factory creation fails because DiagPoint is absent. Fixtures use literal coarse/fine field values, nonzero domain indices and physical origin, and MPI ownership away from the I/O rank.
- [x] Implement parameter validation, finest-cell selection on every prepare, GPU-safe value transfer and MPI collection, covered-cell rejection, and checked ASCII append. Add both build-system entries.
- [x] Run 2D serial, 2D MPI, and 3D MPI tests for values and field ordering, regridding/fallback, boundaries, restart append/deduplication, bad positions/parameters, covered cells, and incompatible output histories. All 60 tests passed after migration, including omission of a volume-fraction field. Test failures use executable exit status and output files, not source text matching.
- [x] Document usable input and exact output semantics in PelePhysics `Docs/sphinx/Utility.rst` and PeleC `Docs/sphinx/InputFiles.rst`. Format and diff checks passed; focused code review found no remaining substantive issues. The actual 2D `PeleC-Sod` build and smoke tests passed with 1 and 2 MPI ranks, covering initialization, selected state/derived fields, two probe intervals, and checkpoint restart append.

**Resolved review finding:** `PeleC::do_diagnostics()` labeled initialization with `parent->cumTime() + dt`, colliding with the first completed step. The production smoke test reproduced the failure before the fix; using the actual state time in `Source/PeleC.cpp` passed initialization, completed-step, and checkpoint-restart checks.

**Validation boundary:** Tests used GCC 11.5 and the current checked-out AMReX/SUNDIALS sources. After transient `/tmp` build caches disappeared, final builds used a persistent WSL task cache; the tested implementation, build inputs, and tests were verified byte-for-byte against the workspace. No GPU execution or real EB geometry run was performed. The EB test supplies synthetic volume fractions through the real sampling path. WarpX FieldProbe was read as a design reference; no WarpX-specific field or particle dependency was imported.
