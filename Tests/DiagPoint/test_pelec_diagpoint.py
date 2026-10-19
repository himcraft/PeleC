"""Run a built 2D PeleC-Sod executable through init, timesteps and restart."""

import argparse
import math
from pathlib import Path
import subprocess
import tempfile


def read_history(path):
    lines = path.read_text().splitlines()
    assert len([line for line in lines if line.startswith("#")]) == 2
    return [[float(value) for value in line.split()]
            for line in lines if line and not line.startswith("#")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    args = parser.parse_args()
    command = [str(args.executable.resolve())]
    if args.ranks > 1:
        command = ["mpiexec", "-n", str(args.ranks)] + command
    with tempfile.TemporaryDirectory(prefix="pelec-diagpoint-smoke-") as directory:
        work = Path(directory)
        (work / "inputs").write_text("""
max_step = 2
stop_time = 1
geometry.is_periodic = 0 0
geometry.coord_sys = 0
geometry.prob_lo = 0 0
geometry.prob_hi = 1 0.25
amr.n_cell = 32 8
amr.max_level = 0
amr.blocking_factor = 4
amr.max_grid_size = 8
amr.plot_files_output = 0
amr.checkpoint_files_output = 1
amr.check_file = chk
amr.check_int = 2
amr.v = 0
eb2.geom_type = all_regular
pelec.lo_bc = Hard SlipWall
pelec.hi_bc = Hard SlipWall
pelec.do_hydro = 1
pelec.do_react = 0
pelec.diffuse_vel = 0
pelec.diffuse_temp = 0
pelec.diffuse_spec = 0
pelec.fixed_dt = 0.001
pelec.init_shrink = 1
pelec.sum_interval = 0
pelec.v = 0
prob.p_l = 1
prob.u_l = 0
prob.rho_l = 1
prob.p_r = 0.1
prob.u_r = 0
prob.rho_r = 0.125
prob.idir = 1
prob.frac = 0.5
pelec.diagnostics = probe1 probe2
pelec.probe1.type = DiagPoint
pelec.probe1.position = 0.125 0.125
pelec.probe1.field_names = pressure Temp density
pelec.probe1.int = 1
pelec.probe1.file = probe1.dat
pelec.probe1.volume_fraction_field = vfrac
pelec.probe2.type = DiagPoint
pelec.probe2.position = 0.875 0.125
pelec.probe2.field_names = pressure density
pelec.probe2.int = 2
pelec.probe2.file = probe2.dat
pelec.probe2.volume_fraction_field = vfrac
""")

        def run(*overrides):
            result = subprocess.run(command + ["inputs", *overrides], cwd=work,
                                    capture_output=True, text=True, timeout=90)
            if result.returncode:
                raise AssertionError(result.stdout + result.stderr)

        run()
        first = read_history(work / "probe1.dat")
        assert [row[0] for row in first] == [0, 1, 2], first
        assert [row[1] for row in first] == [0, 0.001, 0.002], first
        assert math.isclose(first[0][7], 1.0, rel_tol=1e-12), first[0]
        assert first[0][9] == 1.0, first[0]
        assert [row[0] for row in read_history(work / "probe2.dat")] == [0, 2]
        checkpoint = work / "chk00002"
        assert checkpoint.is_dir(), list(work.iterdir())
        run("max_step=4", "amr.restart=chk00002")
        appended = read_history(work / "probe1.dat")
        assert [row[0] for row in appended] == [0, 1, 2, 3, 4], appended
        for step, row in enumerate(appended):
            assert math.isclose(row[1], step * 0.001, abs_tol=1e-15), row
            assert math.isclose(row[7], 1.0, rel_tol=1e-12), row
            assert math.isclose(row[9], 1.0, rel_tol=1e-12), row
        assert [row[0] for row in read_history(work / "probe2.dat")] == [0, 2, 4]
        print(f"PASS: 2D PeleC init, two probes, interval, restart; {args.ranks} rank(s)")


if __name__ == "__main__":
    main()
