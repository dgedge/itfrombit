"""Backend-independent fermion terms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

from .charges import ChargeVector, Mode, add_charge, neg_charge, zero_charge

OpKind = Literal["create", "annihilate", "number"]


@dataclass(frozen=True)
class FermionOp:
    kind: OpKind
    site: int

    def hc(self) -> "FermionOp":
        if self.kind == "create":
            return FermionOp("annihilate", self.site)
        if self.kind == "annihilate":
            return FermionOp("create", self.site)
        return self


@dataclass(frozen=True)
class FermionTerm:
    coeff: complex
    ops: tuple[FermionOp, ...]
    label: str = ""

    def hc(self) -> "FermionTerm":
        return FermionTerm(
            coeff=self.coeff.conjugate(),
            ops=tuple(op.hc() for op in reversed(self.ops)),
            label=f"{self.label}.hc" if self.label else "",
        )


def op_charge(op: FermionOp, mode_by_site: dict[int, Mode]) -> ChargeVector:
    charge = mode_by_site[op.site].charge
    if op.kind == "create":
        return charge
    if op.kind == "annihilate":
        return neg_charge(charge)
    return zero_charge(len(charge))


def term_charge(term: FermionTerm, modes: Iterable[Mode]) -> ChargeVector:
    mode_by_site = {mode.index: mode for mode in modes}
    if not mode_by_site:
        raise ValueError("no modes supplied")
    rank = len(next(iter(mode_by_site.values())).charge)
    total = zero_charge(rank)
    for op in term.ops:
        total = add_charge(total, op_charge(op, mode_by_site))
    return total


def charge_violations(terms: Iterable[FermionTerm], modes: Iterable[Mode]) -> list[FermionTerm]:
    mode_list = list(modes)
    rank = len(mode_list[0].charge)
    zero = zero_charge(rank)
    return [term for term in terms if term_charge(term, mode_list) != zero]


def add_onsite_number(terms: list[FermionTerm], site: int, coeff: complex, label: str) -> None:
    if abs(coeff) > 0:
        terms.append(FermionTerm(coeff, (FermionOp("number", site),), label))


def add_hopping(
    terms: list[FermionTerm],
    i: int,
    j: int,
    coeff: complex,
    label: str,
    *,
    plus_hc: bool = True,
) -> None:
    if i == j:
        raise ValueError("hopping sites must differ")
    if abs(coeff) == 0:
        return
    term = FermionTerm(coeff, (FermionOp("create", i), FermionOp("annihilate", j)), label)
    terms.append(term)
    if plus_hc:
        terms.append(term.hc())
