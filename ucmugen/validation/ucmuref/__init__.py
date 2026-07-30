"""Validation support for UCMuGen, the C++ port of the UCMuon generator.

The Fortran generator in ``src/generator`` is the reference implementation.
This package drives it reproducibly (:mod:`ucmuref.fortran`), compares samples
with null-calibrated tolerances (:mod:`ucmuref.stats`), and defines the matrix
of configurations the port must reproduce (:mod:`ucmuref.cases`).
"""

from .fortran import (
    ANGULAR_NAMES,
    COLUMNS,
    SPECTRUM_NAMES,
    SURFACE_NAMES,
    GenConfig,
    GeneratorError,
    derived,
    parse,
    run,
)

__all__ = [
    "ANGULAR_NAMES", "COLUMNS", "SPECTRUM_NAMES", "SURFACE_NAMES",
    "GenConfig", "GeneratorError", "derived", "parse", "run",
]
