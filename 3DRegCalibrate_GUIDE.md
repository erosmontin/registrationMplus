# 3DRegCalibrate Guide

`3DRegCalibrate` is a metric calibration utility for the registrationMplus
multi-metric cost function. It does not register images or write a transformed
image. Instead, it evaluates the Mplus metric around the identity
`Similarity3DTransform` and reports how large each enabled metric term is over a
controlled set of rotations, translations, and scale changes.

Use it when you want a data-driven starting point for metric weights before
running `3DRegSimilarity`, `3DRegAffine`, `3DRegAffineMultiLevel`, or
`3DRegBsplines`.

## Quick Start

Calibrate MI and NGF value and derivative scales:

```bash
3DRegCalibrate \
  --fixedimage fixed.nii.gz \
  --movingimage moving.nii.gz \
  --metrics "1,1,0,0,0,0" \
  --compute-derivative true \
  --output-csv calibrate_mi_ngf.csv
```

The command prints a summary table and writes one row per perturbation sample to
`calibrate_mi_ngf.csv`.

In `3DRegCalibrate`, `--metrics` is a selector mask, not a weight array. The
order matches the registration tools: `MI,NGF,MSE,GD,NC,NMI[,Label]`. Any
non-zero value means "include this metric in calibration"; the tool internally
evaluates selected metrics with unit weight so the reported ranges are intrinsic
metric scales.

Use `--metric-derivatives` only when the derivative selector should differ from
the value selector:

```bash
3DRegCalibrate \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  --metrics "1,1,0,0,0,0" \
  --metric-derivatives "1,0,0,0,0,0" \
  --compute-derivative true \
  --output-csv calibrate_mi_only_deriv.csv
```

Boolean options use explicit values, for example `--compute-derivative true`
and `--metric-overlap false`.

## What It Does

The tool:

1. Loads the fixed and moving images.
2. Optionally resamples both images to `--workingresolution`.
3. Builds an identity `Similarity3DTransform` centered on the fixed image.
4. Initializes the same `itk::Mplus` composite metric used by the registration
   executables.
5. Sweeps one transform degree of freedom at a time:
   - `rotX`, `rotY`, `rotZ` in degrees.
   - `transX`, `transY`, `transZ` in mm.
   - `scale` as an isotropic scale factor.
6. Records unit-weight sub-metric values, absolute values, derivative norms,
   and derivative absolute ranges at every sample.
7. Prints relative value and derivative weights that can be copied back into
   `--alpha`, `--lambda`, `--rho`, and the corresponding derivative flags.

The default sweep evaluates:

```text
3 rotation axes * 9 samples
3 translation axes * 9 samples
1 scale axis * 5 samples
= 59 metric evaluations
```

## Important Assumptions

`3DRegCalibrate` sweeps around the identity transform. It works best when the
input images are already in approximately the same physical space, or when the
moving image has already been roughly aligned by a prior registration stage.

If the images are far apart, first run a coarse similarity or affine
registration, resample the moving image, and calibrate on that pre-aligned pair.

Only selected metrics are initialized. Prefer `--metrics` for calibration
selection. The individual flags (`--alpha`, `--lambda`, `--rho`, and friends)
are still accepted for compatibility, but in this tool they now behave as
include/exclude selectors; their numeric values are not used as calibration
weights.

The summary uses MI as the reference when MI is selected and has non-zero scale.
If MI is not selected, the first active metric with non-zero scale becomes the
reference. This makes MI-free calibration possible without manually recomputing
ratios from the CSV.

## Interpreting the Summary

The summary includes these columns:

| Column | Meaning |
|--------|---------|
| `Metric` | Sub-metric name: MI, NGF, MSE, GD, NC, NMI, or Label |
| `Value?` | Whether the metric value was selected |
| `Min`, `Max`, `Range` | Unit-weight metric value range |
| `AbsRange` | Range of `abs(metric value)`, used first for value scaling |
| `MeanAbs` | Mean absolute unit-weight metric value |
| `RelValueW` | Suggested value weight relative to the reference metric |
| `Deriv?` | Whether the derivative was selected, shown with `--compute-derivative true` |
| `MeanDerNorm` | Mean per-metric derivative L2 norm |
| `MeanDerAbs` | Mean absolute derivative component value |
| `DerAbsRange` | Mean absolute derivative component range |
| `RelDerivW` | Suggested derivative weight relative to the reference derivative |

Example interpretation:

```text
Metric  AbsRange  RelValueW  MeanDerNorm  MeanDerAbs  RelDerivW
MI      0.12000   1.00000    0.80000      0.18000     1.00000
NGF     0.04000   3.00000    0.20000      0.05000     4.00000
MSE     2.40000   0.05000    8.00000      1.60000     0.10000
```

This means NGF's value changed about one third as much as MI over the same
sweep, so a starting NGF value weight near `3.0` would make its value scale
comparable to MI. Its derivative norm is one quarter of MI, so the derivative
weight could start near `4.0`.

Treat the suggested weights as starting points, not final registration settings.
After calibration, run a registration on representative cases and inspect
alignment quality, convergence, and label/Dice behavior if labels are available.

## CSV Output

Use `--output-csv path.csv` to save all samples.

The CSV columns are:

```text
axis,param_value,
MI,MI_abs,MI_range,MI_abs_range,MI_rel_value_weight,
MI_deriv_norm,MI_deriv_mean_abs,MI_deriv_abs_range,MI_rel_deriv_weight,
...
Total,Total_abs,TotalGradNorm
```

`axis` identifies the sweep axis. `param_value` is degrees for rotation axes,
millimeters for translation axes, and the scale factor for the scale axis.

Each metric is reported as a unit-weight value plus an absolute-value column.
The `_range`, `_abs_range`, `_rel_value_weight`, and
`_rel_deriv_weight` columns repeat the final summary values on every row so the
CSV is self-contained. Derivative columns are `0` unless
`--compute-derivative true` is used.

## Common Recipes

### Multimodal MI + NGF

```bash
3DRegCalibrate \
  -f fixed_ct.nii.gz \
  -m moving_mr.nii.gz \
  --metrics "1,1,0,0,0,0" \
  --metric-sampling "0.1,0.1,0.05,0.05,0.05,0.05" \
  --compute-derivative true \
  --output-csv ct_mr_mi_ngf.csv
```

### Single-Modality MSE + NC with MI Reference

```bash
3DRegCalibrate \
  -f baseline_mr.nii.gz \
  -m followup_mr.nii.gz \
  --metrics "1,0,1,0,1,0" \
  --metric-sampling "0.1,0.05,0.1,0.05,0.1,0.05" \
  --compute-derivative true \
  --output-csv mr_mse_nc.csv
```

### Compare All Image Metrics

```bash
3DRegCalibrate \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  --metrics "1,1,1,1,1,1" \
  --compute-derivative true \
  --normalizemse true \
  --normalizegd true \
  --output-csv all_metrics.csv
```

### Label Metric Calibration

```bash
3DRegCalibrate \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  --fixedlabelmap fixed_labels.nii.gz \
  --movinglabelmap moving_labels.nii.gz \
  --metrics "1,0,0,0,0,0,1" \
  --compute-derivative true \
  --labelsamples 0.1 \
  --labeldistmax 20 \
  --output-csv mi_label.csv
```

For per-label weights, use the legacy map syntax accepted by this executable:

```bash
--labelkappavec "1:1.0,2:0.5,3:0.25"
--labelkappaderivvec "1:1.0,2:0.5,3:0.25"
```

### Faster Coarse Calibration

```bash
3DRegCalibrate \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  --metrics "1,1,0,0,0,0" \
  --workingresolution 4,4,4 \
  --rot-steps 5 \
  --trans-steps 5 \
  --scale-steps 3 \
  --metric-sampling "0.03,0.03,0.03,0.03,0.03,0.03"
```

### Focus on a Region of Interest

```bash
3DRegCalibrate \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  --focusroi brain_mask.nii.gz \
  --metrics "1,1,0,0,0,0" \
  --output-csv roi_calibration.csv
```

## Option Reference

### Required Inputs

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--fixedimage` | `-f` | required | Fixed/reference image |
| `--movingimage` | `-m` | required | Moving image evaluated against the fixed image |

### Calibration Selection Arrays

| Option | Default | Description |
|--------|---------|-------------|
| `--metrics` | empty | Selector mask `MI,NGF,MSE,GD,NC,NMI[,Label]`; non-zero includes the metric value at unit calibration weight |
| `--metric-derivatives` | empty | Selector mask for derivative calibration; defaults to the `--metrics` selector when `--metrics` is provided |
| `--metric-sampling` | empty | Sampling fractions `MI,NGF,MSE,GD,NC,NMI[,Label]`; overrides individual sampling flags |

### Optional Masks and Labels

| Option | Default | Description |
|--------|---------|-------------|
| `--focusroi` | `N` | Fixed-image mask limiting metric evaluation |
| `--fixedlabelmap` | `N` | Fixed label map |
| `--movinglabelmap` | `N` | Moving label map |
| `--labelkappa` | `0.0` | Legacy label value selector when `--metrics` is not used |
| `--labelkappaderiv` | `0.0` | Legacy label derivative selector when `--metric-derivatives` is not used |
| `--labelkappavec` | empty | Per-label value weights as `label:weight,label:weight` |
| `--labelkappaderivvec` | empty | Per-label derivative weights as `label:weight,label:weight` |
| `--labelsamples` | `0.05` | Label samples as a fraction `(0,1]` or absolute count `>1` |
| `--labeldistmax` | `20.0` | Clamp label signed distances in mm |
| `--labelhuber` | `false` | Use the registration-style robust Huber label loss |
| `--labelhuberdelta` | `0.25` | Huber delta for normalized label residuals |
| `--labelnarrowband` | `false` | Evaluate label loss only near either label boundary |
| `--labelbandwidth` | `5.0` | Narrow-band half-width in mm |

### Legacy Individual Selectors

| Option | Default | Description |
|--------|---------|-------------|
| `--alpha` | `1.0` | MI value selector |
| `--alphaderivative` | `1.0` | MI derivative selector |
| `--lambda` | `0.0` | NGF value selector |
| `--lambdaderivative` | `0.0` | NGF derivative selector |
| `--nu` | `0.0` | MSE value selector |
| `--nuderivative` | `0.0` | MSE derivative selector |
| `--rho` | `0.0` | GD value selector |
| `--rhoderivative` | `0.0` | GD derivative selector |
| `--yota` | `0.0` | NC value selector |
| `--yotaderivative` | `0.0` | NC derivative selector |
| `--sigma` | `0.0` | NMI value selector |
| `--sigmaderivative` | `0.0` | NMI derivative selector |

### Metric Sampling and Parameters

| Option | Short | Default | Description |
|--------|-------|---------|-------------|
| `--mattesnumberofbins` | `-b` | `64` | MI histogram bins |
| `--mattespercentage` | `-p` | `0.05` | MI sampling fraction |
| `--ngfpercentage` | | `0.05` | NGF sampling fraction |
| `--msepercentage` | | `0.05` | MSE sampling fraction |
| `--gdpercentage` | | `0.05` | GD sampling fraction |
| `--ncpercentage` | | `0.05` | NC sampling fraction |
| `--nmibins` | | `64` | NMI histogram bins |
| `--nmipercentage` | | `0.05` | NMI sampling fraction |
| `--normalizemse` | | `false` | Normalize MSE |
| `--normalizegd` | | `false` | Normalize GD |
| `--NGFevaluator` | | `0` | NGF evaluator (`0=scalar`) |
| `--ngfprecompute` | | `false` | Registration-compatible NGF gradient precompute option |
| `--etavaluefixed` | | `-1` | Fixed-image NGF eta (`-1=auto`) |
| `--etavaluemoving` | | `-1` | Moving-image NGF eta (`-1=auto`) |
| `--ngfspacing` | | `4,4,4` | NGF finite-difference spacing in mm |

### Sweep Configuration

| Option | Default | Description |
|--------|---------|-------------|
| `--rot-max` | `20.0` | Sweep rotations from `-rot-max` to `+rot-max` degrees |
| `--rot-steps` | `9` | Rotation samples per axis; odd values include zero |
| `--trans-max` | `10.0` | Sweep translations from `-trans-max` to `+trans-max` mm |
| `--trans-steps` | `9` | Translation samples per axis; odd values include zero |
| `--scale-max` | `0.1` | Sweep scale from `1-scale-max` to `1+scale-max` |
| `--scale-steps` | `5` | Scale samples; set `0` to skip scale |
| `--compute-derivative` | `false` | Also call `GetValueAndDerivative` and report gradient norms |

### Runtime and Output

| Option | Default | Description |
|--------|---------|-------------|
| `--numberofthreads` | `2` | CPU threads used by sub-metrics |
| `--bsplinecaching` | `false` | B-spline weight caching; `false` is recommended for calibration |
| `--workingresolution` | `0,0,0` | Resample fixed/moving images to this spacing; `0,0,0` keeps native spacing |
| `--output-csv` | empty | Optional CSV output path |
| `--metric-overlap` | `true` | Restrict metric evaluation to the fixed/moving overlap |
| `--help` | | Print command-line help |

## Applying the Relative Weights

The summary reports separate weights for values and derivatives. A typical
workflow is:

1. Enable the metrics you want to compare with `--metrics`.
2. Run `3DRegCalibrate --compute-derivative true`.
3. Use `RelValueW` as starting `--alpha`, `--lambda`, `--nu`, `--rho`,
   `--yota`, `--sigma`, or `--labelkappa`.
4. Use `RelDerivW` as the corresponding derivative weight.
5. Validate on representative cases.

Example:

```bash
# Suppose calibration suggests NGF RelValueW=2.5 and RelDerivW=4.0.
3DRegAffine \
  -f fixed.nii.gz \
  -m moving.nii.gz \
  -o registered.nii.gz \
  --alpha 1.0 \
  --alphaderivative 1.0 \
  --lambda 2.5 \
  --lambdaderivative 4.0
```

## Troubleshooting

| Symptom | Likely cause | What to try |
|---------|--------------|-------------|
| All non-MI metrics are zero | They were not selected | Use `--metrics`, for example `"1,1,0,0,0,0"` |
| Derivative norms are zero | Derivative evaluation is disabled or derivatives were not selected | Use `--compute-derivative true` and select derivatives with `--metrics` or `--metric-derivatives` |
| Calibration is slow | Too many samples or high native resolution | Use fewer sweep steps, lower sampling fractions, or `--workingresolution` |
| Relative weights look extreme | The metric range is tiny or the images are not roughly aligned | Pre-align the moving image and rerun calibration |
| Label metric is zero | Label maps are missing, label was not selected, or labels do not overlap enough | Provide both label maps and include the 7th `--metrics` value |
| Command rejects a boolean option | Boost expects an explicit bool value | Use `true` or `false` after the option |
