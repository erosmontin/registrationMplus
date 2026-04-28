#!/bin/bash
# B-spline deformable registration (step 2 after rigid)
# Uses rigid .tfm as initialisation via --transformin
# Metrics/derivatives/sampling passed as arrays

LIST=/data/ARTICLE/bodymodels/selected_ixi_multisite_30.txt
T1DIR=/data/MYDATA/IXI/IXI-T1
OPT=/data/ARTICLE/bodymodels/reconnew/recT1_1_5
M=/data/ARTICLE/bodymodels/reconnew/1.5/Duke_tse_ETL1_T1w_2mm/volume.nii.gz
NT=$(nproc)

ROUT=$OPT/rigid
BOUT=$OPT/bspline
mkdir -p "$BOUT"

while IFS= read -r G; do
    [ -n "$G" ] || continue

    F="$T1DIR/$G"
    PT=${G%.nii.gz}

    # --- outputs ---
    O="$BOUT/$G"
    OVD="$BOUT/$PT.mha"
    OTFM="$BOUT/$PT.tfm"
    OSNAP="$OPT/bspline_snap/$PT"

    # --- rigid transform from step 1 ---
    RIGID_TFM="$ROUT/$PT.tfm"

    mkdir -p "$OSNAP"

    # check input image exists
    if [ ! -f "$F" ]; then
        echo "missing input: $F"
        continue
    fi

    # skip if bspline output already done
    if [ -f "$O" ]; then
        echo "skipped (already done): $O"
        continue
    fi

    # require rigid transform from step 1
    if [ ! -f "$RIGID_TFM" ]; then
        echo "skipped (no rigid tfm): $RIGID_TFM"
        continue
    fi

    echo "registering $PT ..."

    3DRegBsplines \
      -f "$F" \
      -m "$M" \
      -o "$O" \
      -v "$OVD" \
      -T "$OTFM" \
      -W "$RIGID_TFM" \
      --numberofthreads "$NT" \
      \
      --metrics          "1.0,0.5,0.1,0,0,0" \
      --metric-derivatives "1.0,0,0,0,0,0" \
      --metric-sampling  "0.05,0.05,0.05,0,0,0" \
      \
      --gridresolution 20 \
      --overlappadding 3 \
      -I 1000 \
      -F 1e7 \
      --normalizemse 1 \
      --normalizegd 1 \
      --ngfprecompute 1 \
      --ngfspacing "4,4,4" \
      \
      --snapshotdir "$OSNAP" \
      --snapshotevery 50 \
      --metricoverlap 1 \
      --verbose 1

    if [ $? -eq 0 ]; then
        echo "done: $O"
    else
        echo "FAILED: $PT" >&2
    fi

done < "$LIST"
