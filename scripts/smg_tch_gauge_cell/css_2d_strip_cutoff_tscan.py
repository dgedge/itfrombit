#!/usr/bin/env python3
"""
Larger-cutoff 2D-strip mirror-SMG gap + t-response scan.

This is the explicit next step named in css_2d_strip_mirror_fock_lanczos.py's own verdict:
  "the next strip run must use a LARGER LOCAL CUTOFF and verify a VISIBLE t-RESPONSE."

At states_per_charge=2 (the minimal cutoff that contains a mirror excitation) the full gap is
positive but the hopping t does not visibly move it -- so a positive gap there cannot be
distinguished from a truncation artifact. This driver pushes the cutoff up (3, 4, ...) and
scans the hopping t at each cutoff. The Gauss-law-projected dimension grows fast with the
cutoff (local_dim = 3 * states_per_charge per vertex, 6 vertices), so this is a big-RAM
sparse-Lanczos job -- run it on a many-core / high-memory box.

READ:
  * gap STAYS POSITIVE as the cutoff grows           -> gap is stable, not a small-cutoff fluke
  * gap VISIBLY MOVES with t (not t-inert)           -> a real DYNAMICAL gap, not a kinematic
                                                        truncation artifact
  both together  -> evidence the mirror SMG gap is real on the (still truncated, candidate-map)
                    spatial strip.  t-inert OR collapsing-with-cutoff -> artifact.

HONEST SCOPE (unchanged from the parent script): this is a finite-volume, charge-resolved,
Gauss-law-truncated strip with a SUBSTITUTED corner->codeword map (items 77/97/102 open). A
clean result is evidence for spatial gap stability, NOT the 3+1D constructive theorem.

Usage (on the box):
  cd <repo>/python_code
  CUTOFFS=2,3,4  python css_2d_strip_cutoff_tscan.py      # or edit CUTOFFS/TS/BETA below
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import css_2d_strip_mirror_fock_lanczos as strip

BETA = float(os.environ.get("BETA", "1.0"))
TS = [float(x) for x in os.environ.get("TS", "0.0,0.1,0.2,0.4,0.8").split(",")]
CUTOFFS = [int(x) for x in os.environ.get("CUTOFFS", "2,3,4").split(",")]

print(f"beta={BETA}  cutoffs={CUTOFFS}  t-scan={TS}", flush=True)
print("(gap should stay >0 across cutoffs AND respond to t for a real dynamical gap)", flush=True)
for spc in CUTOFFS:
    print(f"\n### states_per_charge = {spc} ###", flush=True)
    print(f"  {'t':>5} {'dim':>11} {'nnz':>13} {'full_gap':>10} {'matter_gap':>11}  build+solve", flush=True)
    gaps = []
    for t in TS:
        t0 = time.time()
        try:
            m = strip.StripHamiltonian(beta=BETA, hopping=t, states_per_charge=spc)
            herm = strip.hermiticity_probe(m)
            fg, mg = strip.strip_matter_gap(m)
            gaps.append((t, fg, mg))
            print(f"  {t:>5g} {m.dim:>11d} {len(m.data):>13d} {fg:>10.6g} "
                  f"{strip.plaquette.format_gap(mg):>11}  ({time.time()-t0:.0f}s, herm {herm:.0e})", flush=True)
        except MemoryError:
            print(f"  {t:>5g}  MemoryError -- states_per_charge={spc} too large for this box's RAM", flush=True)
            break
        except Exception as exc:
            print(f"  {t:>5g}  FAILED: {type(exc).__name__}: {exc}", flush=True)
            break
    if len(gaps) >= 2:
        fgs = [g[1] for g in gaps]
        spread = max(fgs) - min(fgs)
        print(f"  -> full-gap range over t: [{min(fgs):.4g}, {max(fgs):.4g}]  (t-spread {spread:.4g}; "
              f"{'RESPONDS to t' if spread > 1e-3 else 'T-INERT (artifact-like)'})", flush=True)

print("\nVERDICT KEY:", flush=True)
print("  gap stays >0 as cutoff grows AND t-spread grows/visible  -> mirror SMG gap looks REAL", flush=True)
print("  gap collapses with cutoff OR stays t-inert at large cutoff -> truncation artifact", flush=True)
print("DONE", flush=True)
