"""
    :file:     modal.py
    :author:   Zhu Dengda (zhudengda@mail.iggcas.ac.cn)
    :contributors: Yuan Yusong (yuanyusong25@gmail.com)
    :date:     2026-06-30

    Python端使用的面波频散相关类。
"""

from __future__ import annotations

import numpy as np

from .c_structures import NPCT_REAL_TYPE

__all__ = [
    "PyDispersion",
]


class PyDispersion:
    def __init__(self, wave_type, freqs, c, ciref, cnum):
        wave_type = _normalize_wave_type(wave_type)
        freqs = np.asarray(freqs, dtype=NPCT_REAL_TYPE)
        c = np.asarray(c, dtype=NPCT_REAL_TYPE)
        ciref = np.asarray(ciref, dtype=np.int64)
        cnum = np.asarray(cnum, dtype=np.int64)

        if freqs.ndim != 1:
            raise ValueError("freqs must be a 1-D array.")
        if c.ndim != 2:
            raise ValueError("c must be a 2-D array.")
        if ciref.shape != c.shape:
            raise ValueError("ciref must have the same shape as c.")
        if c.shape[0] != freqs.size:
            raise ValueError("c and freqs have inconsistent lengths.")
        if cnum.shape != (freqs.size,):
            raise ValueError("cnum must have the same length as freqs.")

        self.wave_type = wave_type
        self.freqs = freqs.copy()
        self.c = c.copy()
        self.ciref = ciref.copy()
        self.cnum = cnum.copy()
        self.modes = np.arange(self.c.shape[1], dtype=np.int64)

    def mode(self, n):
        if n < 0 or n >= self.c.shape[1]:
            raise IndexError(f"mode {n} is out of range.")
        return self.freqs.copy(), self.c[:, n].copy()


def _normalize_wave_type(wave_type):
    wave_type = str(wave_type).strip().lower()
    if wave_type in ("rayleigh", "rayl", "r"):
        return "Rayleigh"
    if wave_type in ("love", "l"):
        return "Love"
    raise ValueError("wave_type must be 'Rayleigh' or 'Love'.")
