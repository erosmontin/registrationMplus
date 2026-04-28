# Snapshot Resolution Control Guide

## Overview
Snapshot resolution control allows you to customize how iteration snapshots are resampled to 2D slices for visualization during registration. Three command-line flags control this behavior.

In addition to resolution/interpolation controls, snapshots now enforce the same transform-space path used by registration so panels are always rendered in the correct fixed-image frame.

## Settings

### 1. `--snapshotspacing` (Recommended — Absolute Spacing)
Sets the **absolute pixel spacing in millimeters (mm)** for the resampled snapshot slice.

**Default:** `0.0` (disabled; uses scale-based method instead)

**Usage:**
```bash
3DRegBsplines --snapshotspacing 0.5 ...
3DRegBsplines --snapshotspacing 1.0 ...
3DRegBsplines --snapshotspacing 2.0 ...
```

**Example — To set snapshot to 0.5 mm per pixel:**
```bash
3DRegBsplines \
  --snapshotspacing 0.5 \
  --snapshotinterp 1 \
  ...
```

**When to use:**
- You want consistent, predictable snapshot pixel size regardless of working resolution
- You're comparing snapshots across registrations with different working resolutions
- You want high-quality snapshots at sub-millimeter precision

**Note:** When `--snapshotspacing` > 0, it **overrides** `--snapshotscale` completely.

---

### 2. `--snapshotscale` (Relative Scaling)
Scales the **minimum working voxel spacing** by a factor.

**Default:** `1.0` (use minimum working spacing as-is)

**Formula:**
```
snapshot_spacing = min_working_voxel_spacing × snapshotscale
```

**Usage:**
```bash
3DRegBsplines --snapshotscale 0.5 ...   # Half as large (2× zoom)
3DRegBsplines --snapshotscale 1.0 ...   # Same size (default)
3DRegBsplines --snapshotscale 2.0 ...   # Twice as large (2× shrink)
```

**Example:**
If working resolution is `[2mm, 2mm, 2mm]` (min = 2mm):
- `--snapshotscale 0.5` → snapshot spacing = 1mm/pixel
- `--snapshotscale 1.0` → snapshot spacing = 2mm/pixel
- `--snapshotscale 2.0` → snapshot spacing = 4mm/pixel

**When to use:**
- You want snapshots to adapt to the working resolution automatically
- Processing multiple datasets with different native resolutions

**When NOT to use:**
- If you need absolute, predictable snapshot sizes → use `--snapshotspacing` instead

---

### 3. `--snapshotinterp` (Interpolation Quality)
Controls resampling interpolator quality.

**Default:** `1` (cubic B-spline)

**Options:**
- `0` — Linear interpolation (faster, lower quality)
- `1` — Cubic B-spline interpolation (recommended; default)

**Usage:**
```bash
3DRegBsplines --snapshotinterp 0 ...   # Linear (fast)
3DRegBsplines --snapshotinterp 1 ...   # Cubic (smooth, default)
```

**When to use:**
- Linear (`0`): Fast prototyping, large batches where speed matters more than smoothness
- Cubic (`1`): Production use, publication-quality snapshots, visual assessment

---

## Decision Matrix

| Scenario | Settings |
|----------|----------|
| **Absolute mm/pixel, smooth** | `--snapshotspacing 0.5 --snapshotinterp 1` |
| **Same as above, fast** | `--snapshotspacing 0.5 --snapshotinterp 0` |
| **Adaptive to working res, smooth** | `--snapshotscale 0.5 --snapshotinterp 1` |
| **Default (no change from current)** | (omit all flags) |
| **Largest snapshots, smooth** | `--snapshotscale 2.0 --snapshotinterp 1` |

---

## Space Consistency Rules

### Working-space vs original-space behavior
- If snapshot target spacing is similar to working spacing, snapshots are sourced from working images.
- If snapshot target spacing is finer (for example, working 2.0 mm and snapshot 0.5 mm), snapshots are sourced from original images and resampled directly to snapshot spacing.

### Transform composition (critical)
- Snapshots apply the same transform chain as registration for moving-image rendering.
- In `3DRegBsplines`, if `--transformin` is a linear transform (Affine/Rigid/Similarity), registration pre-warps moving data before B-spline optimization.
- Snapshot rendering now composes the transforms in the same order:
  - moving->prewarped via initial linear transform
  - prewarped->fixed via current B-spline transform
- Result: snapshot moving panel, checkerboard, and grid overlays stay in the same space as the optimizer.

---

## Code Locations

### Binaries that support these flags:
- `3DRegAffine`
- `3DRegAffineMultiLevel`
- `3DRegSimilarity`
- `3DRegBsplines`

### Implementation files:
- **CLI parsing:** `src/3DRegistration/3DRegBsplines/src/3DRegBsplines.cxx` (line ~200)
- **Observer:** `src/includes/registrationUtils.h` — `IterationSnapshotObserver` class
- **Resampling:** `src/includes/registrationUtils.h:ResampleSliceIsotropic()` method

### Related configuration:
- `examples/configs/config_schema.json` — JSON schema documentation

---

## Technical Details

### Snapshot Creation Pipeline
1. **Build display frame** — Use fixed-image geometry; if requested spacing is finer than working spacing, create a finer snapshot grid in fixed space
2. **Resample moving with transform chain** — Apply current registration transform, and if present, the initial pre-warp transform from `--transformin`
3. **Resample fixed source** — Match the same display frame so panels are aligned
4. **Extract slice + intensity mapping** — Convert to display-ready slices
5. **Composite** — Overlay fixed, moving, difference/grid panels
4. **Save** — Write to disk (NIfTI or PNG format)

### Default Behavior (no flags)
- Uses **minimum voxel spacing** from working resolution
- Applies **cubic B-spline** interpolation
- No extra scaling applied

### With `--snapshotspacing 0.5`
- All snapshots have **exactly 0.5 mm/pixel**
- Smooth cubic resampling applied
- Independent of working resolution

---

## Examples

### Example 1: Production snapshots at 0.5 mm resolution
```bash
./3DRegBsplines \
  --fixed fixed.nii \
  --moving moving.nii \
  --transform bspline \
  --metric mplus \
  --snapshotspacing 0.5 \
  --snapshotinterp 1 \
  --output result.nii \
  --outputtransform transform.h5
```

**Result:** High-quality, consistent 0.5 mm/pixel snapshots regardless of working resolution.

### Example 2: Fast prototyping with adaptive scaling
```bash
./3DRegBsplines \
  --fixed fixed.nii \
  --moving moving.nii \
  --transform bspline \
  --metric mplus \
  --snapshotscale 0.5 \
  --snapshotinterp 0 \
  --output result.nii
```

**Result:** Snapshots scaled 2× (zoomed) relative to working resolution, using fast linear interpolation.

### Example 3: Default (no snapshot customization)
```bash
./3DRegBsplines \
  --fixed fixed.nii \
  --moving moving.nii \
  --output result.nii
```

**Result:** Snapshots use minimum working voxel spacing, smooth cubic interpolation.

---

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| Snapshots are too small (pixelated) | Working resolution is coarse (large mm/voxel) | Use `--snapshotspacing` with finer value (e.g., 0.5) |
| Snapshots are too large (huge files) | `--snapshotspacing` too small | Increase to 1.0 or 2.0 |
| Snapshots blur/smooth differently | Interpolator changed | Use `--snapshotinterp 1` for consistency |
| Snapshot moving/fixed panels look in different spaces | Initial transform path not matched | Ensure latest build: snapshot renderer now composes `--transformin` pre-warp + current transform |
| Performance slow during registration | Snapshots being computed every iteration | Reduce frequency in code or disable (future feature) |

---

## Version History

- **Apr 2026:** Snapshot space-consistency update
  - Snapshots can source from original images at fine target spacing
  - Snapshot moving render composes initial linear pre-warp (`--transformin`) and current registration transform
  - Fixed-space alignment maintained for moving/fixed/checker/grid panels
- **Apr 2026:** Snapshot resolution control implemented with three flags
  - `--snapshotscale` (relative scaling via factor)
  - `--snapshotspacing` (absolute mm/pixel)
  - `--snapshotinterp` (linear vs cubic interpolation)
