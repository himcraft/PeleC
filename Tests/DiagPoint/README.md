# PeleC integration tests for DiagPoint

The shared implementation and standalone AMReX tests live in PelePhysics:
`Submodules/PelePhysics/Source/Utility/Diagnostics/DiagPoint.{H,cpp}` and
`Submodules/PelePhysics/Testing/Diagnostics/DiagPoint`.

This directory provides a CMake entry point for those tests:

```sh
cmake -S Tests/DiagPoint -B tmp_build_dir/diagpoint-2d \
  -DAMReX_SPACEDIM=2 -DAMReX_MPI=ON -DAMReX_OMP=OFF
cmake --build tmp_build_dir/diagpoint-2d -j 8
ctest --test-dir tmp_build_dir/diagpoint-2d --output-on-failure
```

To check production integration, build the 2D `PeleC-Sod` target, then run:

```sh
python3 Tests/DiagPoint/test_pelec_diagpoint.py /path/to/PeleC-Sod
python3 Tests/DiagPoint/test_pelec_diagpoint.py /path/to/PeleC-Sod --ranks 2
```

The smoke test checks initialization, state/derived fields, two probe intervals,
and append after checkpoint restart in a temporary directory.
