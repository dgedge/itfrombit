"""TeNPy entry point.

The local environment used to build this package does not have TeNPy installed.
The backend is therefore written against the documented TeNPy 1.x model API but
is not locally runtime-validated here. Production DMRG should be run in an
environment with TeNPy, SciPy, NumPy, and h5py installed.
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any

import numpy as np

from .model_3450 import (
    FLAVOR_CHARGES,
    FLAVOR_NAMES,
    Model3450Spec,
    build_3450_terms,
    initial_product_state_3450,
)
from .models import ModelSpec, build_model_terms, mode_index
from .terms import FermionOp, FermionTerm, charge_violations


def require_tenpy():
    try:
        import tenpy  # type: ignore
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "TeNPy is not installed. Install TeNPy in the production environment "
            "before running DMRG. The exact-ED guards can run without it."
        ) from exc
    return tenpy


def op_name(op: FermionOp) -> str:
    if op.kind == "create":
        return "Cd"
    if op.kind == "annihilate":
        return "C"
    if op.kind == "number":
        return "N"
    raise ValueError(f"unknown op kind: {op.kind}")


def tenpy_lattice_index(site: int, index_mode: str) -> list[int]:
    if index_mode == "flat":
        return [site, 0]
    if index_mode == "3450":
        return [site // 4, site % 4]
    raise ValueError(f"unknown TeNPy site index mode: {index_mode!r}")


def tenpy_local_term(term: FermionTerm, *, index_mode: str = "flat") -> list[tuple[str, list[int]]]:
    """Convert a backend term to TeNPy's add_local_term convention."""

    return [(op_name(op), tenpy_lattice_index(op.site, index_mode)) for op in term.ops]


def spec_from_model_params(model_params: Any) -> ModelSpec:
    raw = model_params.get("smg_spec", None)
    if raw is None:
        raw = {
            key: model_params.get(key, getattr(ModelSpec(), key))
            for key in asdict(ModelSpec()).keys()
        }
    return ModelSpec(**raw)


def normalize_conserve(value: str | None) -> str | None:
    if value is None:
        return None
    if str(value).lower() in ("none", "false", "0", ""):
        return None
    return value


def build_tenpy_model_class():
    require_tenpy()
    from tenpy.models.model import CouplingMPOModel  # type: ignore
    from tenpy.networks.site import FermionSite  # type: ignore

    class SMGTermListModel(CouplingMPOModel):
        """Finite spinless-fermion chain generated from `ModelSpec`."""

        def init_sites(self, model_params):
            conserve = normalize_conserve(model_params.get("conserve", "N"))
            return FermionSite(conserve=conserve)

        def init_terms(self, model_params):
            spec = spec_from_model_params(model_params)
            modes, terms = build_model_terms(spec)
            violations = charge_violations(terms, modes)
            if violations:
                labels = ", ".join(term.label for term in violations[:5])
                raise ValueError(f"charge-conservation gate failed: {labels}")
            for term in terms:
                self.add_local_term(
                    term.coeff,
                    tenpy_local_term(term),
                    category=term.label or "smg_term",
                )

    return SMGTermListModel


def build_generic_tenpy_model_class():
    require_tenpy()
    from tenpy.models.model import CouplingMPOModel  # type: ignore
    from tenpy.networks.site import FermionSite  # type: ignore

    class GenericTermListModel(CouplingMPOModel):
        """Finite spinless-fermion chain from an explicit term payload."""

        def init_sites(self, model_params):
            conserve = normalize_conserve(model_params.get("conserve", None))
            return FermionSite(conserve=conserve)

        def init_terms(self, model_params):
            terms = model_params["fermion_terms"]
            index_mode = model_params.get("site_index_mode", "flat")
            for term in terms:
                self.add_local_term(
                    term.coeff,
                    tenpy_local_term(term, index_mode=index_mode),
                    category=term.label or "fermion_term",
                )

    return GenericTermListModel


def build_tenpy_model(spec: ModelSpec, *, bc_mps: str = "finite"):
    model_cls = build_tenpy_model_class()
    params = {
        "L": spec.total_modes,
        "bc_MPS": bc_mps,
        "bc_x": "open",
        "conserve": "N",
        "smg_spec": asdict(spec),
    }
    return model_cls(params)


def build_tenpy_term_model(
    total_sites: int,
    terms: list[FermionTerm],
    *,
    conserve: str | None = None,
    bc_mps: str = "finite",
    lattice: Any | None = None,
    site_index_mode: str = "flat",
):
    model_cls = build_generic_tenpy_model_class()
    params = {
        "fermion_terms": terms,
        "site_index_mode": site_index_mode,
    }
    if lattice is None:
        params["L"] = total_sites
        params["bc_MPS"] = bc_mps
        params["bc_x"] = "open"
        params["conserve"] = conserve
    else:
        if lattice.N_sites != total_sites:
            raise ValueError(f"lattice has {lattice.N_sites} sites, expected {total_sites}")
        params["lattice"] = lattice
    return model_cls(params)


def build_charged_fermion_site(q1: int, q2: int):
    """A spinless fermion site conserving the reference Q1,Q2 and parity."""

    require_tenpy()
    from tenpy.linalg import np_conserved as npc  # type: ignore
    from tenpy.networks.site import Site  # type: ignore

    chinfo = npc.ChargeInfo([1, 1, 2], ["Q1", "Q2", "parity"])
    leg = npc.LegCharge.from_qflat(chinfo, [[0, 0, 0], [q1, q2, 1]])
    jw = np.array([[1.0, 0.0], [0.0, -1.0]])
    c = np.array([[0.0, 1.0], [0.0, 0.0]])
    cd = np.array([[0.0, 0.0], [1.0, 0.0]])
    n = np.array([[0.0, 0.0], [0.0, 1.0]])
    site = Site(leg, ["empty", "full"], sort_charge=True, JW=jw, C=c, Cd=cd, N=n)
    site.need_JW_string |= {"C", "Cd", "JW"}
    site.hc_ops["C"] = "Cd"
    site.hc_ops["Cd"] = "C"
    site.hc_ops["JW"] = "JW"
    site.hc_ops["N"] = "N"
    site.charge_to_JW_parity = np.array([0, 0, 1], dtype=int)
    return site


def build_3450_charged_lattice(spec: Model3450Spec):
    """Build a 4-flavor unit-cell lattice with the paper's conserved charges."""

    require_tenpy()
    from tenpy.models.lattice import Lattice  # type: ignore

    unit_cell = [
        build_charged_fermion_site(*FLAVOR_CHARGES[flavor])
        for flavor in FLAVOR_NAMES
    ]
    return Lattice([spec.lattice_orbitals], unit_cell, bc="open", bc_MPS="finite")


def is_3450_charge_conserve(conserve: str | None) -> bool:
    return str(conserve).lower() in {"q3450", "3450", "q1q2", "charges"}


def initial_product_state(spec: ModelSpec) -> list[str]:
    """Half-filled simple seed: physical modes filled, mirror modes empty."""

    state = ["empty" for _ in range(spec.total_modes)]
    for cell in range(spec.cells):
        for flavor in range(spec.flavors):
            state[mode_index(spec, cell, "physical", flavor)] = "full"
    return state


def run_ground_dmrg(
    spec: ModelSpec,
    *,
    chi_max: int = 500,
    max_sweeps: int = 20,
    svd_min: float = 1e-10,
    mixer: bool = True,
) -> dict[str, Any]:
    """Run finite DMRG for the generated term-list model."""

    require_tenpy()
    from tenpy.algorithms import dmrg  # type: ignore
    from tenpy.networks.mps import MPS  # type: ignore

    model = build_tenpy_model(spec)
    psi = MPS.from_product_state(
        model.lat.mps_sites(),
        initial_product_state(spec),
        bc=model.lat.bc_MPS,
    )
    dmrg_params = {
        "mixer": mixer,
        "max_sweeps": max_sweeps,
        "trunc_params": {"chi_max": chi_max, "svd_min": svd_min},
    }
    info = dmrg.run(psi, model, dmrg_params)
    energy = info.get("E", None) if isinstance(info, dict) else None
    return {
        "energy": energy,
        "dmrg_info": info,
        "psi": psi,
        "model": model,
        "mirror_bilinear_m2": mirror_bilinear_m2_tenpy(psi, spec),
    }


def run_ground_dmrg_terms(
    total_sites: int,
    terms: list[FermionTerm],
    initial_state: list[str],
    *,
    chi_max: int = 500,
    max_sweeps: int = 20,
    svd_min: float = 1e-10,
    mixer: bool = True,
    conserve: str | None = None,
    max_trunc_err: float | None = None,
    combine: bool = True,
    lattice: Any | None = None,
    site_index_mode: str = "flat",
) -> dict[str, Any]:
    """Run finite DMRG for an explicit term-list model.

    The 3-4-5-0 interaction changes particle number, so the default TeNPy site
    cannot conserve `N` for that model. Multi-U(1) TeNPy charges should be added
    later for performance; this no-conservation path is the correctness-first
    runnable backend.
    """

    require_tenpy()
    from tenpy.algorithms import dmrg  # type: ignore
    from tenpy.networks.mps import MPS  # type: ignore

    model = build_tenpy_term_model(
        total_sites,
        terms,
        conserve=conserve,
        lattice=lattice,
        site_index_mode=site_index_mode,
    )
    mps_kwargs: dict[str, Any] = {}
    if lattice is not None:
        mps_kwargs["unit_cell_width"] = len(model.lat.unit_cell)
    psi = MPS.from_product_state(
        model.lat.mps_sites(),
        initial_state,
        bc=model.lat.bc_MPS,
        **mps_kwargs,
    )
    dmrg_params = {
        "mixer": mixer,
        "max_sweeps": max_sweeps,
        "max_trunc_err": max_trunc_err,
        "combine": combine,
        "trunc_params": {"chi_max": chi_max, "svd_min": svd_min},
    }
    info = dmrg.run(psi, model, dmrg_params)
    energy = info.get("E", None) if isinstance(info, dict) else None
    return {"energy": energy, "dmrg_info": info, "psi": psi, "model": model}


def run_ground_dmrg_3450(
    spec: Model3450Spec,
    *,
    chi_max: int = 500,
    max_sweeps: int = 20,
    svd_min: float = 1e-10,
    mixer: bool = True,
    conserve: str | None = None,
    max_trunc_err: float | None = None,
    combine: bool = True,
) -> dict[str, Any]:
    modes, terms = build_3450_terms(spec)
    violations = charge_violations(terms, modes)
    if violations:
        labels = ", ".join(term.label for term in violations[:5])
        raise ValueError(f"3-4-5-0 charge-conservation gate failed: {labels}")
    use_3450_charges = is_3450_charge_conserve(conserve)
    lattice = build_3450_charged_lattice(spec) if use_3450_charges else None
    return run_ground_dmrg_terms(
        spec.total_modes,
        terms,
        initial_product_state_3450(spec),
        chi_max=chi_max,
        max_sweeps=max_sweeps,
        svd_min=svd_min,
        mixer=mixer,
        conserve=None if use_3450_charges else conserve,
        max_trunc_err=max_trunc_err,
        combine=combine,
        lattice=lattice,
        site_index_mode="3450" if use_3450_charges else "flat",
    )


def pair_term_dagger_pair(
    left_pair: tuple[int, int],
    right_pair: tuple[int, int],
) -> list[tuple[str, int]]:
    i, j = left_pair
    k, l = right_pair
    return [("Cd", i), ("Cd", j), ("C", l), ("C", k)]


def mirror_bilinear_m2_tenpy(psi: Any, spec: ModelSpec) -> float:
    """Finite-size pair-bilinear order proxy from the MPS.

    This deliberately mirrors the ED diagnostic:

        m2 = L^-2 sum_xy <O_x^dag O_y>,  O_x=sum_{a<b} c_{x,a} c_{x,b}.

    TeNPy's `expectation_value_term(..., autoJW=True)` is used for the
    fermionic strings.
    """

    pairs_by_cell: list[list[tuple[int, int]]] = []
    for cell in range(spec.cells):
        pairs = []
        for flavor_a in range(spec.flavors):
            for flavor_b in range(flavor_a + 1, spec.flavors):
                pairs.append(
                    (
                        mode_index(spec, cell, "mirror", flavor_a),
                        mode_index(spec, cell, "mirror", flavor_b),
                    )
                )
        pairs_by_cell.append(pairs)

    total = 0.0 + 0.0j
    for left_pairs in pairs_by_cell:
        for right_pairs in pairs_by_cell:
            for left_pair in left_pairs:
                for right_pair in right_pairs:
                    term = pair_term_dagger_pair(left_pair, right_pair)
                    total += psi.expectation_value_term(term, autoJW=True)
    return float((total / (spec.cells * spec.cells)).real)
