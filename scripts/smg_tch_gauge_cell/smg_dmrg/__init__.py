"""1+1D SMG DMRG prototype package.

The package exposes backend-independent term construction plus small exact
diagonalisation guards. TeNPy is optional and imported lazily by the CLI.
"""

from .charges import ChargeVector, Mode, build_flavor_cartan_modes
from .model_3450 import Model3450Spec, build_3450_terms, free_spectrum_3450
from .models import ModelSpec, build_model_terms
from .terms import FermionOp, FermionTerm

__all__ = [
    "ChargeVector",
    "FermionOp",
    "FermionTerm",
    "Mode",
    "Model3450Spec",
    "ModelSpec",
    "build_3450_terms",
    "build_flavor_cartan_modes",
    "build_model_terms",
    "free_spectrum_3450",
]
