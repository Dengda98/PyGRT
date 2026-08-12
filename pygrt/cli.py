"""
    :file:     cli.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :date:     2026-8-11

    调用 grt 命令行程序的工具

"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
from typing import Iterable, Optional, Sequence, Union


def find_grt() -> str:
    """Return the nearest available ``grt`` executable path."""
    bundled = Path(__file__).resolve().parent / "C_extension" / "bin" / "grt"
    if bundled.is_file():
        return str(bundled)

    from_path = shutil.which("grt")
    if from_path:
        return from_path

    raise FileNotFoundError(
        "The grt executable was not found. Build it with `make -C pygrt/C_extension CC=gcc-16 -j`."
    )


def run_grt(
    args: Sequence[object],
    *,
    cwd: Optional[Union[str, os.PathLike]] = None,
    print_log: bool = True,
) -> None:
    """
    Run one ``grt`` module and raise a detailed error on failure.

    :param    args:      Arguments passed to the ``grt`` command.
    :param    cwd:       Working directory used to run the command.
    :param    print_log: If true, stream ``grt`` stdout/stderr to the terminal.
                         If false, capture them and attach any output to the
                         raised error on failure.
    """
    command = [find_grt(), *(str(arg) for arg in args)]
    if print_log:
        completed = subprocess.run(command, cwd=cwd, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"grt command failed with exit code {completed.returncode}: {' '.join(command)}")
        return

    completed = subprocess.run(command, cwd=cwd, check=False, text=True, capture_output=True)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        message = f"grt command failed with exit code {completed.returncode}: {' '.join(command)}"
        if detail:
            message += f"\n{detail}"
        raise RuntimeError(message)


def format_float(value: float) -> str:
    """
    Format a numerical CLI value without unnecessary trailing zeros.

    :param    value: Numerical value to format.

    :return: Formatted command-line value.
    """
    return format(float(value), ".15g")


def format_range(values: Iterable[float], name: str) -> str:
    """
    Format three values for a ``-X`` or ``-Y`` option.

    :param    values: Three values corresponding to one CLI coordinate option.
    :param    name:   Name used in the validation error message.

    :return: Slash-separated command-line value.
    """
    try:
        values = list(values)
    except TypeError:
        raise ValueError(f"{name} must contain exactly three values.")
    if len(values) != 3:
        raise ValueError(f"{name} must contain exactly three values.")

    return "/".join(format_float(value) for value in values)
