"""Named bounds for the pulse-chase suite -- committed BEFORE the runs they gate.

The suite's acceptance rule (mirrored from `canon.PARITY_TOL`): a threshold chosen
after seeing the number is a rationalization, not a gate. Every "meaningful bound" a
discrimination or robustness check needs lives here as a named constant with its
rationale, so no check carries an inline literal whose value was reverse-engineered
from a passing run.

None of these are canonical-basis values (those live in `canon.py` and may never be
transcribed). These are behavioural bounds on ECSPr itself, so they belong to the
suite, not to the basis -- but they follow the same commit-before-you-look discipline.

Solver/atom tolerances are imported from the engine, never restated:
  * ecspr_solver.SELFTEST_TOL   (1e-9)   -- star/production two-terminal R_eff
  * ecspr_solver.REFF_EPS/IEFF_EPS (1e-12) -- star/production conductance floor
  * ecspr_atom_graph.ATOM_REFF_EPS (1e-12) -- ATOM-lane conductance floor (distinct)
  * canon.EPS (1e-12), canon.PARITY_TOL (1e-9), canon.GPU_CPU_TOL (1e-6)
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# Discrimination bounds (I4, II7) -- "meaningful", not merely > 1.
# ---------------------------------------------------------------------------
# The star manufactures conductance between co-participants that transfer no atom
# to each other (cofactor theft, input->input leak). The magnitude is set by the
# stolen partner's TOTAL atom count, not the transferred count: on MNXR106432 the
# zero-carbon pyruvate-NADPH channel runs ~21x the real one-carbon pyruvate-CO2
# channel. A "meaningful" artifact must clear a factor that float noise or a near-1
# coincidence cannot reach; 2x is the floor we require the star to exceed and the
# atom graph to erase (atom graph puts the two atoms in different components -> 0).
STAR_LEAK_RATIO_MIN = 2.0

# ---------------------------------------------------------------------------
# permute_pairs control (II8) -- the shuffle must collapse the effect.
# ---------------------------------------------------------------------------
# Shuffling which-atom-goes-where within each reaction (preserving counts and
# mappability) must destroy most of the measured effect, or the lane is reporting
# mappability rather than chemistry. We require the shuffled mean |Ieff effect| to
# fall below half the unshuffled -- a real chemical signal loses more than that when
# its atom identities are scrambled.
PERMUTE_COLLAPSE_FRAC = 0.5

# ---------------------------------------------------------------------------
# Far-distance robustness (III1, III2).
# ---------------------------------------------------------------------------
# How many of the longest testable axes per element to exercise as far edges. Pinned
# so "far coverage" is a stated count, not whatever happened to be reachable.
FAR_EDGE_COUNT = 5
# A far-edge base R_eff above this is too close to the 1e-12 conductance clamp for the
# long-path Ieff to be signal rather than float noise. 1/1e6 = 1e-6 sits six orders
# above the clamp, so anything below this bound is safely resolved.
FAR_REFF_MAX = 1.0e6

# ---------------------------------------------------------------------------
# Strict-inequality margin for toy property checks (II2, II5, II9, I4).
# ---------------------------------------------------------------------------
# Toy Ieff values are O(1); a strict ">" claim ("parallel routes exceed either
# alone", "dropping a labelled atom lowers Ieff") must clear this absolute margin to
# count, so float wobble at the 1e-15 level never decides a red/green verdict.
STRICT_EPS = 1.0e-9

# ---------------------------------------------------------------------------
# Directed-conductance throttling (IV1) -- direction must actually bite.
# ---------------------------------------------------------------------------
# An irreversible reaction is a diode: it conducts forward (g+) and is throttled
# backward (g- = ratio * g+). For a single irreversible step the two-terminal
# conductance ratio C_fwd / C_rev is ~1/ratio. The directionality ensemble clamps a
# strongly-favoured reaction to ratio ~ 1e-2 or below, so a real throttle clears ~100x;
# we require the directed solver to show at least this factor between forward and reverse
# on a toy irreversible chain. Below it, the "directed" solve is not actually directional
# -- float noise or a near-symmetric ratio cannot manufacture a 10x gap. The companion
# claim (IV2) is the opposite bound: at ratio == 1 the gap must VANISH (parity), so the
# two checks pin the directed solver from both sides.
DIRECTED_THROTTLE_MIN = 10.0

# The backward/forward ratio of the irreversible reaction in IV1's toy chain. NOT a basis
# value -- a fixture knob, chosen small enough that the series throttle (~1/ratio) clears
# DIRECTED_THROTTLE_MIN with room (~50x here). Named rather than inlined so it is not a bare
# magic literal, and deliberately kept off the CLEAN-floor literal that the transcribed-
# numbers tripwire guards (a different, unrelated basis constant).
IV1_TOY_IRREV_RATIO = 0.02
