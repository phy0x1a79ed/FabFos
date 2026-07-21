# Sockeye — groundwork, and what is left (T5, NOT done)

**Status: not done.** The container — the part the plan called "the real work" — is
built and pin-verified. The cluster facts below are probed, not assumed. What has
not happened is the run: no ECSPr step has executed on Sockeye through the engine.

Stopping here was the plan's own recommendation if T5 exceeded one chunk: ship the
rest and say so plainly rather than slide into a cluster-plumbing project. It did,
so this says so plainly. Everything below is verified, so the next session starts
from facts instead of rediscovering them.

## Verified (probed 2026-07-14)

| Fact | Value |
|---|---|
| Access | `ssh sockeye` already connected as `txyliu`; no human needed |
| 2FA | workspace `2fa` tooling enrolled (`cwl`, `alliance`) — check it before concluding a human is required |
| Apptainer | `apptainer/1.3.1`, and it is **not** visible until `gcc/7.5.0` is loaded **first, in its own separate `module load`** |
| GPU | partitions `gpu` and `interactive_gpu`, both `gpu:v100:4`, 184320 MB |
| Accounts | `st-shallam-1` (CPU) and `st-shallam-1-gpu` (GPU) — they differ by resource class, so a GPU job charged to the CPU account will not run |
| Scratch | `/scratch/st-shallam-1` |
| Container | `container_builds/main/ecspr`, tagged and pin-asserted **in-image at build time** |

The image's pins are checked by the image itself during `docker build`, so a pin
regression fails the build rather than a queued GPU job six hours later. Both pins
matter, but not for the reason usually given — see that recipe's `load/env.yml`.

## What is left

1. **Get the image onto the cluster.** It is local-only (~13 GB docker). Either push
   to quay and let Apptainer pull, or `dev.sh --sif` and transfer. Note the recorded
   bulk-transfer figure for this link is ~11 MB/s (the often-repeated 5–10 KB/s
   number is stale), so a transfer is minutes, not hours.
2. **Point `containers::ecspr.oci` at whatever the cluster can resolve.** It
   currently holds a registry URI that is not pushed yet.
3. **Deploy the agent** over an `SshSource` rooted on scratch. The library's one
   existing SSH driver is the template.
4. **Run the chain** with `--runtime APPTAINER`, requesting the 32 GB GPU constraint
   and charging the GPU account.

## The thing that is already true

`device`/`dtype` are staged, content-hashed inputs (`ecspr::compute_profile`), so the
GPU path is reachable **through the driver** for the first time — which is precisely
what closing T6 means. It was previously unreachable: `context.params` is populated
only from the runtime resources line, so `context.params.get("device")` was always
its default, which is why the last end-to-end run bypassed the engine with raw
sbatch. `canon.COMPUTE_GPU` stages the GPU profile; a GPU run and a CPU run of the
same step produce different cache keys, so they cannot collide.

The spine already runs this chain in this container locally, and reproduces the
incumbent (see `parity/check_spine_output.py`). Sockeye changes where it runs, not
what it computes.

## Before you run it

`canon.GPU_CPU_TOL` is **already committed**, deliberately, before any GPU run
exists to compare against. GPU float64 and CPU float64 should agree closely but not
bit-exactly. Do not adjust that number after seeing the first result — a threshold
chosen after seeing the number is a rationalization, not a gate. If the measured
difference exceeds it, that is a finding, and it needs a human.

Solve and ablation are CPU-intractable on the dense maps, so a CPU fallback is not a
fallback.
