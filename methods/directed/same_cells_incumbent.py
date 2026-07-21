"""Emit the incumbent (undirected) effect on EXACTLY the cells the directed run scored,
so directed-vs-undirected is a clean same-cells, same-conditions comparison. effect = the
incumbent panel's s = log(G_pert/G_base)."""
import sys
from pathlib import Path
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import benchmark_bridge as bb   # reuses adapt(), FACET_MAP, PANEL_OUT

facet = "A_iECDH10B"
directed = pd.read_csv(HERE / "bench_out" / f"observations_directed_{facet}.tsv", sep="\t")
keys = set(zip(directed.condition_id, directed.edge_id))

frames = []
for arm in ("lof", "gof"):
    df = pd.read_parquet(bb.PANEL_OUT / f"panel_{arm}.parquet")
    df = df[(df.facet == facet) & df.usable & ~df.anchor_incident].copy()
    df["effect"] = df["s"].astype(float)
    frames.append(df)
full = pd.concat(frames, ignore_index=True)
obs = bb.adapt(full, facet)
# restrict to the directed run's exact cells
obs = obs[[(c, e) in keys for c, e in zip(obs.condition_id, obs.edge_id)]]
out = HERE / "bench_out" / f"observations_incumbent_samecells_{facet}.tsv"
obs.to_csv(out, sep="\t", index=False)
print(f"wrote {len(obs):,} incumbent same-cells observations (directed had {len(directed):,}) -> {out}")
