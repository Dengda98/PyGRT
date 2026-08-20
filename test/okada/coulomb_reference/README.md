# Coulomb finite-fault reference data

The five complete Coulomb input files under `input/` are copied from the
Coulomb example directory.  Each result directory contains the corresponding
human-readable Coulomb output files and `coulomb.mat`.

The CI comparison uses `coulomb.mat` as the numerical reference.  Its `U` and
`D` arrays were produced by the direct Coulomb `Okada_DC3D`/`Okada_DC3D0`
calculation and retain full precision.  `Displacement.cou` and
`Displacement_derivatives.cou` are retained for inspection.  The ordinary
Coulomb displacement text output is not used as the strict reference because
its writer rounds to eight decimal places and its point-source GUI output has
the known scale issue.

The reference data do not require MATLAB during CI.  Regenerate them only when
the Coulomb calculation or the selected example inputs are intentionally
updated.
