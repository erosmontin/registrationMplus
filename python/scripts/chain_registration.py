#!/usr/bin/env python3
"""Optimised pyramidal chain registration
==========================================

Similarity → AffineMultiLevel → Bspline (coarse-to-fine)

Speed optimisation checklist (applied automatically by this script)
-------------------------------------------------------------------
1.  **Staged pipeline** — rigid 7-DOF ➜ affine 12-DOF ➜ B-spline deformable.
    Each stage is warm-started (`-W`) from the previous transform so the
    optimizer starts near the solution and wastes no iterations.

2.  **Multi-resolution pyramid** — `3DRegAffineMultiLevel` uses ITK's
    `MultiResolutionPyramidImageFilter`. At each level the images are
    Gaussian-smoothed (σ ∝ 1/shrink) and downsampled ×2. Coarsest level runs
    first ➜ broad convergence basin, fast evaluation.

3.  **Coarse-to-fine B-spline** — `3DRegBsplines` is single-resolution, so
    the script calls it multiple times with decreasing `--gridresolution`
    (e.g. 80 → 50 → 30 mm). Each run warm-starts from the previous transform.

4.  **Random metric sampling** — `--mattespercentage`, `--ngfpercentage`,
    `--msepercentage`, `--ncpercentage`, `--labelsamples` all default to 0.1
    (10 %). For the coarse stages where only global alignment matters, even
    lower (5 %) is fine. Use `--sampling` to set all at once.

5.  **NGF precompute** — `--ngfprecompute 1` is always on. The moving-image
    NGF gradient is computed once and resampled each iteration, instead of
    being re-derived from scratch every evaluation.

6.  **ITK global threading** — `ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS` is set
    from `--threads`, capping **all** ITK filter threads, not just the
    optimizer's.

7.  **No snapshots by default** — snapshot I/O (PNG/NIfTI per iteration) is
    off unless you explicitly add `--snapshotdir … --snapshotevery N` via
    `--extra-*-opts`.

8.  **Convergence factor** — B-spline LBFGS-B uses `-F` (cost-function
    convergence factor). Tighter = slower but more accurate; `1e4` is a good
    balance.

Presets
-------
  --preset fast      Low sampling (5 %), fewer iters, 2 affine levels,
                     coarse B-spline only (80,50). Good for previewing.
  --preset balanced  (default) 10 % sampling, 3 affine levels, 3 B-spline
                     levels (80,50,30).
  --preset accurate  15 % sampling, 4 affine levels, 4 B-spline levels
                     (100,60,35,20), more iterations.

Examples
--------
  # Default balanced run
  python chain_registration.py fixed.nii moving.nii \\
      --out-prefix /g/recit/P04 --bin-dir ./build-v2/bin

  # Fast preview with labels
  python chain_registration.py fixed.nii moving.nii --preset fast \\
      --fixed-label seg_f.nii.gz --moving-label seg_m.nii.gz

  # Accurate with CUDA binaries
  python chain_registration.py fixed.nii moving.nii --preset accurate \\
      --bin-dir ./build-cuda/bin --threads 8

  # Dry-run: print commands without executing
  python chain_registration.py fixed.nii moving.nii --dry-run
"""
from __future__ import annotations

import argparse
import math
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path


# ═══════════════════════════════════════════════════════════════════════════
# Presets
# ═══════════════════════════════════════════════════════════════════════════

PRESETS: dict[str, dict] = {
    "fast": dict(
        sampling=0.05,
        sim_iters=100,  sim_maxstep=4.0,  sim_minstep=0.001,
        aff_iters=150,  aff_levels=2,
        bs_grid_schedule="80,50",  bs_iters=150,  bs_convergence=1e5,
    ),
    "balanced": dict(
        sampling=0.10,
        sim_iters=200,  sim_maxstep=2.0,  sim_minstep=0.0001,
        aff_iters=250,  aff_levels=3,
        bs_grid_schedule="80,50,30",  bs_iters=300,  bs_convergence=1e4,
    ),
    "accurate": dict(
        sampling=0.15,
        sim_iters=350,  sim_maxstep=2.0,  sim_minstep=1e-5,
        aff_iters=400,  aff_levels=4,
        bs_grid_schedule="100,60,35,20",  bs_iters=500,  bs_convergence=1e3,
    ),
}


# ═══════════════════════════════════════════════════════════════════════════
# Helpers
# ═══════════════════════════════════════════════════════════════════════════

def run_cmd(cmd: str, env: dict | None = None, dry: bool = False) -> float:
    """Run *cmd* via subprocess. Returns wall-clock seconds."""
    if dry:
        print(f"  [dry-run] {cmd}")
        return 0.0
    print(f"\n>>> {cmd}\n")
    t0 = time.monotonic()
    subprocess.run(shlex.split(cmd), check=True, env=env)
    return time.monotonic() - t0


def _sampling_opts(args) -> str:
    """Return sampling-percentage flags common to every executable."""
    s = args.sampling
    return (
        f" --mattespercentage {s}"
        f" --ngfpercentage {s}"
        f" --msepercentage {s}"
        f" --ncpercentage {s}"
        f" --labelsamples {s}"
    )


def _label_opts(args) -> str:
    if args.fixed_label and args.moving_label:
        return (
            f" --fixedlabelmap {args.fixed_label}"
            f" --movinglabelmap {args.moving_label}"
            f" --labelkappa {args.label_kappa}"
            f" --labelkappaderiv {args.label_kappa_deriv}"
            f" --labelreport {args.label_report}"
        )
    return ""


def _metric_weight_opts(args) -> str:
    """Return the metric weight flags if explicitly set, else use modality."""
    if args.modality != "custom":
        return f" --modality {args.modality}"
    return (
        f" -a {args.mi_weight} -A {args.mi_weight_deriv}"
        f" -l {args.ngf_weight} -L {args.ngf_weight_deriv}"
        f" -n {args.mse_weight} -N {args.mse_weight_deriv}"
    )


def _common_opts(args) -> str:
    """Flags shared by all stages."""
    return (
        f" --numberofthreads {args.threads}"
        f" --ngfprecompute 1"
        f" --ngfspacing {args.ngf_spacing}"
        f" --derivativemode {args.derivative_mode}"
        + _sampling_opts(args)
        + _label_opts(args)
        + _metric_weight_opts(args)
    )


# ═══════════════════════════════════════════════════════════════════════════
# Stage builders
# ═══════════════════════════════════════════════════════════════════════════

def build_similarity_cmd(bin_dir: Path, args) -> tuple[str, str]:
    out_img = args.out_prefix + "_01_similarity.nii.gz"
    out_tf  = args.out_prefix + "_01_similarity.tfm"
    out_vf  = args.out_prefix + "_01_similarity_vf.nii.gz"
    cmd = (
        f"{bin_dir}/3DRegSimilarity"
        f" -f {args.fixed} -m {args.moving}"
        f" -o {out_img} -v {out_vf} -T {out_tf}"
        f" -I {args.sim_iters}"
        f" -S {args.sim_step} -X {args.sim_maxstep} -R {args.sim_minstep}"
        + _common_opts(args)
    )
    if args.extra_similarity_opts:
        cmd += " " + args.extra_similarity_opts
    return cmd, out_tf


def build_affine_cmd(bin_dir: Path, args, in_tf: str) -> tuple[str, str]:
    out_img = args.out_prefix + "_02_affine.nii.gz"
    out_tf  = args.out_prefix + "_02_affine.tfm"
    out_vf  = args.out_prefix + "_02_affine_vf.nii.gz"
    cmd = (
        f"{bin_dir}/3DRegAffineMultiLevel"
        f" -f {args.fixed} -m {args.moving}"
        f" -o {out_img} -v {out_vf} -T {out_tf}"
        f" -W {in_tf}"
        f" -I {args.aff_iters}"
        f" -U {args.aff_levels}"
        + _common_opts(args)
    )
    if args.extra_affine_opts:
        cmd += " " + args.extra_affine_opts
    return cmd, out_tf


def build_bspline_cmds(
    bin_dir: Path, args, in_tf: str
) -> list[tuple[str, str, float]]:
    """One (cmd, out_tf, grid_mm) per coarse-to-fine level."""
    schedule = [float(g.strip()) for g in str(args.bs_grid_schedule).split(",")]
    steps: list[tuple[str, str, float]] = []
    current_tf = in_tf

    for i, grid_mm in enumerate(schedule):
        tag = f"_03_bspline_L{i}_g{int(grid_mm)}"
        out_img = args.out_prefix + tag + ".nii.gz"
        out_tf  = args.out_prefix + tag + ".tfm"
        out_vf  = args.out_prefix + tag + "_vf.nii.gz"
        cmd = (
            f"{bin_dir}/3DRegBsplines"
            f" -f {args.fixed} -m {args.moving}"
            f" -o {out_img} -v {out_vf} -T {out_tf}"
            f" -W {current_tf}"
            f" --gridresolution {grid_mm}"
            f" -I {args.bs_iters}"
            f" -b {args.bs_blocksize}"
            f" -F {args.bs_convergence}"
            f" -E {args.bs_evaluations}"
            f" -C {args.bs_corrections}"
            + _common_opts(args)
        )
        if args.extra_bspline_opts:
            cmd += " " + args.extra_bspline_opts
        steps.append((cmd, out_tf, grid_mm))
        current_tf = out_tf

    return steps


# ═══════════════════════════════════════════════════════════════════════════
# Profile / dry-run summary
# ═══════════════════════════════════════════════════════════════════════════

def _estimate_bspline_params(image_shape: tuple[int, ...], grid_mm: float,
                             spacing_mm: float = 1.0) -> int:
    """Rough estimate of B-spline DoF for a given grid spacing."""
    n = 1
    for dim_mm in image_shape:
        nodes = int(dim_mm * spacing_mm / grid_mm) + 3  # +3 for cubic border
        n *= max(nodes, 4)
    return n * 3  # 3D vector per node


def print_profile(args):
    """Print a summary of what will run, without executing anything."""
    schedule = [float(g.strip()) for g in str(args.bs_grid_schedule).split(",")]
    n_bs = len(schedule)
    total_stages = 1 + 1 + n_bs

    print()
    print("╔══════════════════════════════════════════════════════════╗")
    print("║           Registration Pipeline — Dry Run               ║")
    print("╠══════════════════════════════════════════════════════════╣")
    print(f"║  Fixed  : {Path(args.fixed).name:>44s} ║")
    print(f"║  Moving : {Path(args.moving).name:>44s} ║")
    print(f"║  Threads: {args.threads:>3d}  (ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS)    ║")
    print(f"║  Sampling: {args.sampling*100:.0f}%  (Mattes, NGF, MSE, NC, Label)       ║")
    print(f"║  NGF precompute: ON                                    ║")
    print(f"║  Snapshots: OFF (add via --extra-*-opts if needed)     ║")
    print("╠══════════════════════════════════════════════════════════╣")
    print(f"║  Stage 1  Similarity       {args.sim_iters:>4d} iters × 7 DOF          ║")
    print(f"║  Stage 2  AffineMultiLevel {args.aff_iters:>4d} iters × {args.aff_levels} levels × 12 DOF  ║")
    for i, g in enumerate(schedule):
        print(f"║  Stage 3.{i}  Bspline g={g:>3.0f}mm  {args.bs_iters:>4d} iters (LBFGS-B)     ║")
    print(f"║                                                        ║")
    print(f"║  Total sub-stages: {total_stages}                                   ║")
    print("╚══════════════════════════════════════════════════════════╝")


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Optimised pyramidal chain registration: "
            "Similarity → AffineMultiLevel (ITK pyramid) → Bspline (coarse→fine)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("fixed",  help="Fixed image filename")
    parser.add_argument("moving", help="Moving image filename")

    # ── I/O ─────────────────────────────────────────────────────────────────
    g_io = parser.add_argument_group("I/O")
    g_io.add_argument("--fixed-label",  default=None, help="Fixed labelmap")
    g_io.add_argument("--moving-label", default=None, help="Moving labelmap")
    g_io.add_argument("--out-prefix",  default="reg_out", help="Output prefix")
    g_io.add_argument("--bin-dir",     default="./build-v2/bin",
                      help="Directory with 3DReg* binaries")

    # ── Preset / mode ──────────────────────────────────────────────────────
    g_mode = parser.add_argument_group("Preset / mode")
    g_mode.add_argument("--preset", choices=["fast", "balanced", "accurate"],
                        default=None,
                        help="Apply a tuned preset (overrides individual params "
                             "unless you also set them explicitly)")
    g_mode.add_argument("--dry-run", action="store_true",
                        help="Print the command plan without executing")

    # ── Global tuning ──────────────────────────────────────────────────────
    g_g = parser.add_argument_group("Global tuning")
    g_g.add_argument("--threads",        type=int,   default=4,
                     help="Caps ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS + "
                          "--numberofthreads for every binary")
    g_g.add_argument("--derivative-mode",type=int,   default=1,
                     help="0=weighted-sum, 1=normalised (RSGD), 2=adaptive")
    g_g.add_argument("--ngf-spacing",    default="5,5,5",
                     help="NGF gradient spacing x,y,z (voxels)")
    g_g.add_argument("--sampling",       type=float, default=0.1,
                     help="Fraction of voxels sampled for EVERY metric "
                          "(Mattes, NGF, MSE, NC, Label). 0.05–0.15 typical.")
    g_g.add_argument("--modality",       default="custom",
                     choices=["multimodal", "singlemodal", "custom"],
                     help="Apply a modality preset for metric weights")

    # ── Metric weights (only used if --modality custom) ────────────────────
    g_w = parser.add_argument_group("Metric weights (--modality custom)")
    g_w.add_argument("--mi-weight",       type=float, default=1.0)
    g_w.add_argument("--mi-weight-deriv", type=float, default=1.0)
    g_w.add_argument("--ngf-weight",      type=float, default=0.2)
    g_w.add_argument("--ngf-weight-deriv",type=float, default=0.2)
    g_w.add_argument("--mse-weight",      type=float, default=0.0)
    g_w.add_argument("--mse-weight-deriv",type=float, default=0.0)

    # ── Label ──────────────────────────────────────────────────────────────
    g_l = parser.add_argument_group("Label metric (all stages)")
    g_l.add_argument("--label-kappa",       type=float, default=0.3)
    g_l.add_argument("--label-kappa-deriv", type=float, default=0.3)
    g_l.add_argument("--label-report",      type=int,   default=5)

    # ── Stage 1 ────────────────────────────────────────────────────────────
    g_s = parser.add_argument_group("Stage 1 — Similarity (RSGD)")
    g_s.add_argument("--sim-iters",   type=int,   default=200)
    g_s.add_argument("--sim-step",    type=float, default=0.001)
    g_s.add_argument("--sim-maxstep", type=float, default=2.0)
    g_s.add_argument("--sim-minstep", type=float, default=0.0001)
    g_s.add_argument("--extra-similarity-opts", dest="extra_similarity_opts",
                     default="")

    # ── Stage 2 ────────────────────────────────────────────────────────────
    g_a = parser.add_argument_group(
        "Stage 2 — AffineMultiLevel (ITK Gaussian pyramid)")
    g_a.add_argument("--aff-iters",  type=int, default=250,
                     help="Iterations per pyramid level")
    g_a.add_argument("--aff-levels", type=int, default=3,
                     help="Pyramid levels (-U). ITK halves resolution + "
                          "Gaussian-smooths at each level. Coarsest first.")
    g_a.add_argument("--extra-affine-opts", dest="extra_affine_opts", default="")

    # ── Stage 3 ────────────────────────────────────────────────────────────
    g_b = parser.add_argument_group("Stage 3 — Bspline coarse-to-fine")
    g_b.add_argument("--bs-grid-schedule", default="80,50,30",
                     help="Comma-separated grid spacings mm, coarse → fine")
    g_b.add_argument("--bs-iters",       type=int,   default=300)
    g_b.add_argument("--bs-blocksize",   type=int,   default=16,
                     help="Mattes histogram bins (-b)")
    g_b.add_argument("--bs-convergence", type=float, default=1e4,
                     help="LBFGS-B convergence factor (-F)")
    g_b.add_argument("--bs-evaluations", type=int,   default=500,
                     help="Max function evaluations per LBFGS-B run (-E)")
    g_b.add_argument("--bs-corrections", type=int,   default=5,
                     help="LBFGS-B corrections (-C)")
    g_b.add_argument("--extra-bspline-opts", dest="extra_bspline_opts", default="")

    args = parser.parse_args(argv)

    # ── Apply preset (overwrite only defaults, not explicit user args) ─────
    if args.preset:
        p = PRESETS[args.preset]
        raw = vars(args)
        # argparse doesn't expose "was this explicitly set", so we key off
        # whether the value still equals the declared default:
        defaults = {a.dest: a.default
                    for a in parser._actions if hasattr(a, "dest")}
        for k, v in p.items():
            if k in defaults and raw.get(k) == defaults.get(k):
                raw[k] = v

    # ── Validate ───────────────────────────────────────────────────────────
    bin_dir = Path(args.bin_dir).resolve()
    if not bin_dir.exists():
        print(f"ERROR: binary directory not found: {bin_dir}", file=sys.stderr)
        sys.exit(2)
    if bool(args.fixed_label) != bool(args.moving_label):
        print("ERROR: supply both --fixed-label and --moving-label, or neither.",
              file=sys.stderr)
        sys.exit(2)

    env = os.environ.copy()
    env["ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS"] = str(args.threads)

    schedule = [float(g.strip()) for g in str(args.bs_grid_schedule).split(",")]
    n_bs = len(schedule)

    # ── Profile / dry-run header ───────────────────────────────────────────
    print_profile(args)
    if args.dry_run:
        print("\nCommands that would be executed:\n")

    timings: dict[str, float] = {}

    # ── Stage 1: Similarity ────────────────────────────────────────────────
    sim_cmd, sim_tf = build_similarity_cmd(bin_dir, args)
    try:
        t = run_cmd(sim_cmd, env=env, dry=args.dry_run)
        timings["similarity"] = t
        if not args.dry_run:
            print(f"[1/3] Similarity done ({t:.1f}s)  →  {sim_tf}")
    except subprocess.CalledProcessError as e:
        print("ERROR: Similarity failed:", e, file=sys.stderr)
        sys.exit(1)

    # ── Stage 2: AffineMultiLevel ──────────────────────────────────────────
    aff_cmd, aff_tf = build_affine_cmd(bin_dir, args, sim_tf)
    try:
        t = run_cmd(aff_cmd, env=env, dry=args.dry_run)
        timings["affine"] = t
        if not args.dry_run:
            print(f"[2/3] AffineMultiLevel done ({t:.1f}s)  →  {aff_tf}")
    except subprocess.CalledProcessError as e:
        print("ERROR: AffineMultiLevel failed:", e, file=sys.stderr)
        sys.exit(2)

    # ── Stage 3: Bspline multi-scale ───────────────────────────────────────
    bs_steps = build_bspline_cmds(bin_dir, args, aff_tf)
    final_bs_tf = aff_tf
    for i, (bs_cmd, bs_tf, grid_mm) in enumerate(bs_steps):
        try:
            t = run_cmd(bs_cmd, env=env, dry=args.dry_run)
            timings[f"bspline_g{int(grid_mm)}"] = t
            final_bs_tf = bs_tf
            if not args.dry_run:
                print(f"[3/3] Bspline L{i} g={grid_mm:.0f}mm ({t:.1f}s)  →  {bs_tf}")
        except subprocess.CalledProcessError as e:
            print(f"ERROR: Bspline L{i} g={grid_mm}mm failed:", e,
                  file=sys.stderr)
            sys.exit(3)

    # ── Summary ────────────────────────────────────────────────────────────
    if args.dry_run:
        print("\n(dry-run mode — nothing was executed)")
    else:
        total = sum(timings.values())
        print()
        print("═" * 60)
        print("All stages complete.")
        print(f"  Final image     : {args.out_prefix}_03_bspline_L{n_bs-1}_g{int(schedule[-1])}.nii.gz")
        print(f"  Final transform : {final_bs_tf}")
        print(f"  Wall time       : {total:.1f}s")
        for stage, dt in timings.items():
            pct = dt / total * 100 if total > 0 else 0
            print(f"    {stage:20s}  {dt:7.1f}s  ({pct:4.1f}%)")
        print("═" * 60)


if __name__ == "__main__":
    main()
