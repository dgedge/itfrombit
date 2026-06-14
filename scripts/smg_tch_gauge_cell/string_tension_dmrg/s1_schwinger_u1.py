#!/usr/bin/env python3
"""
S1 — U(1) Schwinger Kogut-Susskind DMRG: validate the KS-LGT machinery.

Goal: confirm a DMRG of a 1+1D compact-U(1) lattice gauge theory gives a CONFINING (linear)
static quark-antiquark potential V(r), and extract the string tension sigma = slope. This is the
machinery check for the string-tension program (scope.md); U(1) has no tri-colour 2/3 factor, so
this validates the method, not the (4/3)Lambda prediction (that is S2, Z3).

Formulation (gauge-eliminated, open BC; Banuls-Cichy-Cirac-Jansen): integrate the gauge field out
via Gauss law in temporal gauge. Staggered (Kogut-Susskind) fermions on N sites; the electric
energy becomes the running-charge Coulomb term

    H = -w sum_n (c^dag_n c_{n+1} + h.c.)  +  m sum_n (-1)^n n_n  +  J sum_{n=0}^{N-2} L_n^2,
    L_n = l0 + sum_{k<=n} q_k,   q_k = n_k - (1-(-1)^k)/2 + Q^ext_k,   J = g^2/2.

Expanding sum_n L_n^2 (Q^ext are c-numbers):
  * density-density (Q-independent):  2J (N-1-k) n_j n_k  for j<k
  * onsite (Q-dependent via c_n):     [ m(-1)^k + J( (N-1-k) + 2 sum_{n>=k} c_n ) ] n_k
  * constant (Q-dependent):           J sum_n c_n^2     (tracked, added to E)
where c_n = l0 + sum_{k<=n} (Q^ext_k - (1-(-1)^k)/2).

Static potential: V(r) = E_GS(+1 at i_L, -1 at i_R=i_L+r) - E_GS(no external charge).
Confining  <=>  V(r) linear, slope sigma > 0.  Self-checks: vacuum is charge-neutral half-filling.
"""
import sys, numpy as np, warnings
warnings.filterwarnings("ignore")
from tenpy.networks.site import FermionSite
from tenpy.models.model import CouplingMPOModel
from tenpy.algorithms import dmrg
from tenpy.networks.mps import MPS


def bg(k):
    return (1 - (-1) ** k) / 2.0


class Schwinger(CouplingMPOModel):
    default_lattice = "Chain"
    force_default_lattice = True

    def init_sites(self, mp):
        return FermionSite(conserve="N")

    def init_terms(self, mp):
        N = mp["L"]; w = mp["w"]; m = mp["m"]; J = mp["J"]
        Qext = mp.get("Qext", {}) or {}
        l0 = mp.get("l0", 0.0)
        # c_n on links n=0..N-2
        c = [0.0] * (N - 1)
        run = l0
        for n in range(N - 1):
            run += (Qext.get(n, 0.0) - bg(n))
            c[n] = run
        # suffix[k] = sum_{n=k}^{N-2} c_n
        suffix = [0.0] * N
        s = 0.0
        for k in range(N - 1, -1, -1):
            if k <= N - 2:
                s += c[k]
            suffix[k] = s
        # onsite (non-uniform): mass + electric-linear
        for k in range(N):
            coeff = m * (-1) ** k + J * ((N - 1 - k) + 2.0 * suffix[k])
            self.add_onsite_term(coeff, k, "N")
        # hopping (uniform, JW handled by add_coupling)
        self.add_coupling(-w, 0, "Cd", 0, "C", 1, plus_hc=True)
        # density-density Coulomb (non-uniform, no JW for N N)
        for kk in range(N):
            for jj in range(kk):
                self.add_coupling_term(2.0 * J * (N - 1 - kk), jj, kk, "N", "N")
        self.coulomb_const = J * sum(cn * cn for cn in c)


def ground_energy(N, w, m, J, Qext=None, chi=128, sweeps=24):
    Qext = Qext or {}
    M = Schwinger({"L": N, "w": w, "m": m, "J": J, "Qext": Qext, "bc_MPS": "finite"})
    init = ["full" if (k % 2 == 1) else "empty" for k in range(N)]  # neutral staggered vacuum
    psi = MPS.from_product_state(M.lat.mps_sites(), init, bc="finite")
    info = dmrg.run(psi, M, {"mixer": True, "max_sweeps": sweeps,
                             "trunc_params": {"chi_max": chi, "svd_min": 1e-12}})
    Ntot = float(sum(psi.expectation_value("N")))
    return info["E"] + M.coulomb_const, Ntot, float(max(psi.chi))


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    w, m, J = 1.0, 0.5, 0.5
    print(f"S1 Schwinger U(1) KS-DMRG: N={N}, w={w}, m={m}, J=g^2/2={J}", flush=True)
    E0, Ntot, chi = ground_energy(N, w, m, J, {})
    print(f"vacuum: E0={E0:.5f}  N_tot={Ntot:.3f} (expect {N//2})  chi_used={chi:.0f}", flush=True)
    rs, Vs = [], []
    for r in range(1, min(8, N // 2)):
        iL = (N - r) // 2
        iR = iL + r
        Er, _, _ = ground_energy(N, w, m, J, {iL: +1.0, iR: -1.0})
        V = Er - E0
        rs.append(r); Vs.append(V)
        print(f"  r={r:>2}  iL={iL} iR={iR}  V(r)={V:.5f}", flush=True)
    rs = np.array(rs, float); Vs = np.array(Vs, float)
    A = np.vstack([rs, np.ones(len(rs))]).T
    (sigma, b), res, *_ = np.linalg.lstsq(A, Vs, rcond=None)
    ss_tot = np.sum((Vs - Vs.mean()) ** 2)
    r2 = 1 - (np.sum((Vs - (sigma * rs + b)) ** 2) / ss_tot if ss_tot > 0 else 0)
    print(f"\nlinear fit V(r)=sigma*r+b:  sigma={sigma:.4f}  b={b:.4f}  R^2={r2:.4f}", flush=True)
    print(f"sqrt(sigma)={np.sqrt(max(sigma,0)):.4f} (lattice units); slope ~ J={J} expected if confining",
          flush=True)
    print("VALIDATION:", "PASS (linear, sigma>0)" if (sigma > 0.05 and r2 > 0.95)
          else "CHECK (not cleanly linear/confining)", flush=True)


if __name__ == "__main__":
    main()
