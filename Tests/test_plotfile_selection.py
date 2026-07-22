#!/usr/bin/env python3

"""Check regular and small plotfile variable selection."""

from pathlib import Path


def read_plot_variables(plotfile):
    header = plotfile / "Header"
    if not header.is_file():
        raise AssertionError(f"Missing plotfile header: {header}")

    with header.open(encoding="utf-8") as stream:
        stream.readline()
        number_of_variables = int(stream.readline())
        variables = [stream.readline().strip() for _ in range(number_of_variables)]

    if len(variables) != len(set(variables)):
        raise AssertionError(f"Duplicate variables in {header}: {variables}")

    level_header = plotfile / "Level_0" / "Cell_H"
    if not level_header.is_file() or level_header.stat().st_size == 0:
        raise AssertionError(f"Missing or empty MultiFab header: {level_header}")

    data_files = list((plotfile / "Level_0").glob("Cell_D_*"))
    if not data_files or any(path.stat().st_size == 0 for path in data_files):
        raise AssertionError(f"Missing or empty MultiFab data in {plotfile}")

    return set(variables)


def check_plot_variables(plotfile, expected):
    actual = read_plot_variables(plotfile)
    if actual != expected:
        missing = sorted(expected - actual)
        unexpected = sorted(actual - expected)
        raise AssertionError(
            f"Wrong variables in {plotfile}: missing={missing}, "
            f"unexpected={unexpected}"
        )


def main():
    derived_variables = {
        "pressure",
        "x_velocity",
        "Y(H2)",
        "Y(O2)",
        "Y(H2O)",
        "Y(OH)",
    }

    check_plot_variables(Path("plt00000"), derived_variables)
    check_plot_variables(
        Path("smlplt00000"), derived_variables | {"Temp", "heatRelease"}
    )


if __name__ == "__main__":
    main()
