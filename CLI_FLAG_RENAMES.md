# CLI Flag Renames

This document records all flag renames introduced to improve clarity.
Old names are **removed** — update any scripts or pipelines accordingly.

---

## 3DRegBsplines

| Old flag | New flag | Type | Unit | Meaning |
|---|---|---|---|---|
| `--metricoverlap` | `--metric-overlap` | bool | — | Restrict metric to fixed/moving overlap region |
| `--metricpadding` | `--metric-padding-mm` | double | mm | Shrink the overlap region inward on each side |
| `--meshmarginsize` | `--grid-margin-mm` | double | mm | Expand the B-spline grid domain outward |
| `--overlappadding` | `--grid-border-knots` | uint | knots | Extra B-spline knots outside the image border per side |

### Example — before and after

**Before:**
```bash
3DRegBsplines \
  --metricoverlap 1 \
  --metricpadding 30 \
  --meshmarginsize 10 \
  --overlappadding 5
```

**After:**
```bash
3DRegBsplines \
  --metric-overlap 1 \
  --metric-padding-mm 30 \
  --grid-margin-mm 10 \
  --grid-border-knots 5
```

---

## 3DRegAffine / 3DRegSimilarity / 3DRegAffineMultiLevel

| Old flag | New flag | Type | Unit | Meaning |
|---|---|---|---|---|
| `--metricoverlap` | `--metric-overlap` | bool | — | Restrict metric to fixed/moving overlap region |
| `--overlappadding` | `--metric-padding-voxels` | uint | voxels | Shrink the overlap region inward on each side |

> Note: `--grid-margin-mm` and `--grid-border-knots` do not exist in the affine executables — they are B-spline specific.

### Example — before and after

**Before:**
```bash
3DRegAffine \
  --metricoverlap 1 \
  --overlappadding 20
```

**After:**
```bash
3DRegAffine \
  --metric-overlap 1 \
  --metric-padding-voxels 20
```

---

## Focus ROI (3DRegBsplines only — new feature)

| Flag | Type | Default | Meaning |
|---|---|---|---|
| `--focusroi` | string | `N` | Binary or label mask. Non-zero voxels define the region of interest |
| `--focusboost` | double | `0.8` | Fraction of samples drawn from inside the mask (1.0 = hard mask, 0.0 = uniform) |

### Example

```bash
3DRegBsplines \
  --focusroi femur_mask.nii \
  --focusboost 0.8
```

This draws 80 % of metric samples from inside the mask and 20 % from the rest of the image.

---

## Quick migration checklist for agents

Replace these strings in every call to a registration executable:

```
metricoverlap       →  metric-overlap
metricpadding       →  metric-padding-mm        (3DRegBsplines only)
meshmarginsize      →  grid-margin-mm           (3DRegBsplines only)
overlappadding      →  grid-border-knots        (3DRegBsplines)
overlappadding      →  metric-padding-voxels    (3DRegAffine / Similarity / MultiLevel)
```
