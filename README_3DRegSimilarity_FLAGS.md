# 3DRegSimilarity CLI Flag Reference

This document describes the effect of all CLI flags for `3DRegSimilarity`.
It is based on the current implementation in `src/3DRegistration/3DRegAffine/src/3DRegSimilarity.cxx` and the built binary in `bld2/bin/3DRegSimilarity` (v1.2.1).

## Basic Usage

```bash
3DRegSimilarity -f fixed.nii.gz -m moving.nii.gz -o registered.nii.gz [options]
```

Required inputs are `-f`, `-m`, and `-o`. If any are missing, help text is shown and the program exits.

## Important Parsing Note (Booleans)

Boolean options are implemented with `po::value<bool>()`, not switch-style toggles.
That means they require an explicit value:

```bash
--verbose 1
--snapshotgrid 0
--metricoverlap true
```

Using `--verbose` (without value) will throw a command-line parse error.

## 1) General I/O and Control Flags

| Short | Long | Default | Effect |
|---|---|---:|---|
| `-h` | `--help` | n/a | Print help text and exit. |
|  | `--version` | n/a | Print version and exit. |
| `-f` | `--fixedimage` | required | Path to fixed (target) image. |
| `-m` | `--movingimage` | required | Path to moving image. |
| `-o` | `--outputimage` | required | Path to output resampled/registered image. |
| `-v` | `--vfout` | `N` | Write deformation vector field if not `N`. |
| `-T` | `--transformout` | `N` | Write final transform (`.tfm`) if not `N`. |
| `-W` | `--transformin` | `N` | Initialize from transform file if not `N`; otherwise uses identity + center/moments initialization. |
| `-P` | `--dfltpixelvalue` | `0` | Fill value for out-of-bounds voxels during resampling. |
| `-V` | `--verbose` | `0` | Print all parsed options at startup. |
|  | `--numberofthreads` | `2` | Set threads used by registration and metric. |

## 2) Optimizer Flags (Regular Step Gradient Descent)

| Short | Long | Default | Effect |
|---|---|---:|---|
| `-I` | `--maxnumberofiterations` | `1000` | Maximum optimizer iterations. |
| `-S` | `--minimumsteplength` | `0.1` | Minimum step length before stopping. |
| `-X` | `--maximumsteplength` | `1.0` | Initial/maximum step length. |
| `-R` | `--relaxationfactor` | `0.5` | Step reduction factor. |
| `-G` | `--gradientmagnitudetolerance` | `1e-4` | Stop when gradient magnitude is below this tolerance. |

## 3) Composite Metric Weights (Value Terms)

The total objective combines multiple metrics. Setting a weight to `0` effectively turns that value term off.

| Short | Long | Default | Effect |
|---|---|---:|---|
| `-a` | `--alpha` | `1.0` | MI (Mattes Mutual Information) value weight. |
| `-l` | `--lambda` | `1.0` | NGF value weight. |
| `-n` | `--nu` | `1.0` | MSE value weight. |
| `-y` | `--yota` | `0` | NC (normalized correlation) value weight. |
|  | `--rho` | `0` | GD (gradient difference) value weight. |
|  | `--sigma` | `0` | NMI value weight. |
|  | `--labelkappa` | `0` | Label-map value weight. |

## 4) Composite Metric Weights (Derivative Terms)

These scale each metric's contribution to the gradient.

| Short | Long | Default | Effect |
|---|---|---:|---|
| `-A` | `--alphaderivative` | `1.0` | MI derivative weight. |
| `-L` | `--lambdaderivative` | `0` | NGF derivative weight. |
| `-N` | `--nuderivative` | `1.0` | MSE derivative weight. |
| `-Y` | `--yotaderivative` | `0` | NC derivative weight. |
|  | `--rhoderivative` | `0` | GD derivative weight. |
|  | `--sigmaderivative` | `0` | NMI derivative weight. |
|  | `--labelkappaderiv` | `0` | Label-map derivative weight. |

## 5) Sampling / Histogram / NGF Flags

| Short | Long | Default | Effect |
|---|---|---:|---|
| `-p` | `--mattespercentage` | `0.1` | Fraction of voxels sampled for MI. |
| `-b` | `--mattesnumberofbins` | `64` | MI joint-histogram bins. |
|  | `--explicitPDFderivatives` | `0` | Use explicit PDF derivatives for MI. |
|  | `--ngfpercentage` | `0.1` | Fraction of voxels sampled for NGF. |
|  | `--msepercentage` | `0.1` | Fraction of voxels sampled for MSE. |
|  | `--ncpercentage` | `0.1` | Fraction of voxels sampled for NC. |
|  | `--gdpercentage` | `0.1` | Declared as GD sample fraction. |
|  | `--nmipercentage` | `0.1` | Declared as NMI sample fraction. |
|  | `--nmibins` | `64` | NMI joint-histogram bins. |
| `-r` | `--etavaluefixed` | `-1` | NGF fixed-image eta (`-1` enables auto-estimation logic). |
| `-s` | `--etavaluemoving` | `-1` | NGF moving-image eta (`-1` enables auto-estimation logic). |
|  | `--NGFevaluator` | `0` | NGF evaluator variant (`0..4`, validated). |
|  | `--ngfspacing` | `4,4,4` | NGF spacing per axis; must be exactly 3 comma-separated values. |
|  | `--ngfprecompute` | `0` | Precompute moving-image NGF gradient once (faster, approximate). |
|  | `--metricoverlap` | `1` | Compute metric on overlap region between fixed/moving images. |
|  | `--overlappadding` | `20` | Padding (voxels) used around overlap. |
| `-t` | `--fixedimagethreshold` | `-99999999` | If changed from sentinel, only fixed voxels above threshold are considered. |

## 6) Derivative Merge Strategy

| Long | Default | Effect |
|---|---:|---|
| `--derivativemode` | `0` | `0=consistent`, `1=normalized`, `2=main-metric adaptive`. |
| `--mainmetric` | `0` | Main metric index for mode 2 (`0=MI`, `1=NGF`, `2=MSE`, `3=NC`, `4=Label`, `5=GD`, `6=NMI`). |

## 7) Modality Presets

| Long | Default | Effect |
|---|---:|---|
| `--modality` | `custom` | Preset weight profile: `multimodal`, `singlemodal`, or `custom`. |

Preset behavior:

- `multimodal`: defaults to MI + NGF emphasis (`alpha=1`, `alphaderivative=1`, `lambda=0.5`, `lambdaderivative=0.5`, and sets `nu/yota/rho/sigma` plus derivative counterparts to `0` when not explicitly set).
- `singlemodal`: defaults to MSE + NC emphasis (`nu=1`, `nuderivative=1`, `yota=0.5`, `yotaderivative=0.5`, and sets `alpha/lambda/rho/sigma` plus derivative counterparts to `0` when not explicitly set).
- Explicitly provided CLI weights override preset defaults.
- Any other string causes an error and exits.

## 8) Label-Map Flags

| Long | Default | Effect |
|---|---:|---|
| `--fixedlabelmap` | `N` | Fixed label-map path (`N` disables). |
| `--movinglabelmap` | `N` | Moving label-map path (`N` disables). |
| `--labelkappa` | `0` | Global label value weight. |
| `--labelkappaderiv` | `0` | Global label derivative weight. |
| `--labelkappavec` | `""` | Per-label value weights as `L1:w1,L2:w2,...`. |
| `--labelkappaderivvec` | `""` | Per-label derivative weights as `L1:w1,L2:w2,...`. |
| `--labelsamples` | `0.1` | Label metric sampling parameter. |
| `--labelreport` | `1` | Report per-label Dice every N iterations (`0` disables). |

Notes:

- Label metric/observer logic is activated only when both fixed and moving label maps are provided.
- With `labelkappa=0` and `labelkappaderiv=0`, label maps can still be useful for monitoring Dice via `--labelreport`.

## 9) Snapshot / Iteration Visualization Flags

| Long | Default | Effect |
|---|---:|---|
| `--snapshotdir` | `N` | Directory for iteration snapshots (`N` disables snapshots). |
| `--snapshotevery` | `1` | Save snapshot every N iterations. |
| `--snapshotstack` | `0` | If `1`, save full 3D `.nii.gz`; otherwise save 2D snapshot panels. |
| `--snapshotgrid` | `1` | Overlay warped grid on snapshot panels. |
| `--snapshotgridspacing` | `20` | Grid line spacing in voxels. |

## Known Caveats in Current `3DRegSimilarity`

- `--transformin` only accepts `Similarity3DTransform` files; other transform types now produce a clear error instead of a null-pointer crash.
- Help banner text now correctly says "Similarity Registration".

## Minimal Example

```bash
3DRegSimilarity \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  -o registered.nii.gz \
  --modality multimodal \
  -I 500 -S 0.01 -X 1.0 \
  --numberofthreads 8 \
  --verbose 1
```
