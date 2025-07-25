# Image registration based on mPlus
![GitHub last commit](https://img.shields.io/github/last-commit/erosmontin/registrationMplus)
![GitHub issues](https://img.shields.io/github/issues/erosmontin/registrationMplus)

![GitHub forks](https://img.shields.io/github/forks/erosmontin/registrationMplus)
![GitHub stars](https://img.shields.io/github/stars/erosmontin/registrationMplus)


This project implements a multi-metric registration strategy that combines Mutual Information (MI), Normalized Gradient Field (NGF), and Mean Squared Error (MSE) techniques. Developed using the Insight Segmentation and Registration Toolkit (ITK), this method is specifically designed for applications in pediatric oncology.

Pediatric oncology presents a particularly challenging scenario for image registration. Children’s brains undergo significant anatomical changes as they grow, and these changes are further complicated by treatment-induced deformations—such as those caused by hydrocephalus or surgical interventions—which can result in nonuniform and unpredictable alterations in tissue structure. Traditional registration methods often struggle to accurately align images acquired across extended periods, as they may not adequately account for these rapid and heterogeneous changes. To address these challenges, our multi-metric registration strategy leverages ITK's powerful and flexible framework to integrate MI for robust intensity-based alignment, NGF to capture spatial gradients and edge information, and MSE to fine-tune the alignment. This complementary approach is specifically designed for pediatric oncology, where precise image registration is essential for tracking treatment outcomes and correlating radiotherapy dose with neurocognitive effects over long-term follow-up.


For a detailed description of the method, please refer to our article:  
[A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology](https://link.springer.com/article/10.1007/s11517-019-02109-4)

[Publication list](https://biodimensional.com/)



## Features

- **Multi-metric registration:** Combines MI, NGF, and MSE, NC to optimize registration accuracy.
- **Non-rigid 3D registration:** Efficiently handles deformations in brain images.
- **Multithreading support:** Accelerates computation for large 3D datasets.
- **Open source:** Freely available for research and development.

## Installation

### Prerequisites

Make sure your system has the following dependencies installed. You can install them using the commands below:

```bash

sudo apt-get update
sudo apt-get install -y cmake build-essential libinsighttoolkit4-dev
sudo apt-get install -y libpng-dev libjpeg-dev libtiff-dev libdcmtk-dev libfltk1.3-dev libeigen3-dev
sudo apt-get install -y libboost-all-dev
```

## Installation
```bash
git clone https://github.com/erosmontin/registrationMplus.git
cd registrationMplus
mkdir build && cd build
cmake ..
make -j4
sudo make install
```

## Usage
Usage Instructions for Registration Executables

This project provides two registration executables that use the itkMplus metric, each tailored for a different type of registration.

------------------------------------------------------------
1. 3D similarity Registration (itkMplus)
------------------------------------------------------------

This executable performs affine registration using a similarity transform. It leverages the itkMplus metric and a Regular Step Gradient Descent optimizer.
At a minimum, you must supply the fixed image, moving image, output image, and the deformation field output filename.
Additional parameters allow you to fine-tune MI weights, NGF parameters, MSE parameters, and optimizer settings.

Required Options:
  - -f, --fixedimage       : Fixed image filename.
  - -m, --movingimage      : Moving image filename.
  - -o, --outputimage      : Output registered image filename.
  - -v, --vfout            : Deformation field output filename.

Additional Options:
  - --numberofthreads      : Number of threads (default: 2).
  - -a, --alpha            : MI weight (default: 1.0).
  - -A, --alphaderivative   : MI derivative weight (default: 1.0).
  - -l, --lambda           : Lambda value for NGF.
  - -L, --lambdaderivative  : Lambda derivative for NGF.
  - -n, --nu               : Nu value for MSE.
  - -N, --nuderivative      : Nu derivative for MSE.
  - -p, --mattespercentage : Mattes percentage (default: 0.1).
  - -b, --mattesnumberofbins: Number of bins for Mattes (default: 64).
  - Additional options control transform initialization, thresholding, and verbosity.

Example Command:
  3DRegSimilarity \
    -f fixedImage.nii.gz \
    -m movingImage.nii.gz \
    -o registeredImage.nii.gz \
    -v deformationField.nii.gz \
    --numberofthreads 4 \
    -a 1.0 \
    -A 1.0 \
    -l 1.0 \
    -L 0 \
    -n 1.0 \
    -N 1.0 \
    -p 0.1 \
    -b 64

Use the "--help" flag with the executable to see the complete list of options.


------------------------------------------------------------
2. Deformable Registration (B-Spline Registration)
------------------------------------------------------------

This executable (named DeformableRegistration6) performs non-rigid registration using a B‑Spline transform,
optimized with the LBFGSB algorithm. In addition to the standard image inputs/outputs, it accepts options
for configuring the B‑Spline grid (mesh resolution, domain, etc.) and optimizer settings such as convergence factors,
gradient tolerances, and evaluation limits.

Required Options:
  - -f, --fixedimage       : Fixed image filename.
  - -m, --movingimage      : Moving image filename.
  - -o, --outputimage      : Output registered image filename.
  - -v, --vfout            : Deformation field output filename.

Additional Options:
  - --numberofthreads      : Number of threads (default: 2).
  - --gridresolution       : Mesh resolution in mm (e.g., 50).
  - --bsplinecaching       : Enable caching of B‑Spline weights (default: true).
  - -a, --alpha            : MI weight.
  - -A, --alphaderivative   : MI derivative.
  - -l, --lambda           : Lambda value for NGF.
  - -L, --lambdaderivative  : Lambda derivative for NGF.
  - -n, --nu               : Nu value for MSE.
  - -N, --nuderivative      : Nu derivative for MSE.
  - -p, --mattespercentage : Mattes percentage.
  - -b, --mattesnumberofbins: Number of bins for Mattes.
  - Optimizer-specific options include:
      --costfunctionconvergencefactor : Cost function convergence factor.
      --projectedgradienttolerance      : Projected gradient tolerance.
      --numberofevaluations             : Maximum number of evaluations.
      --numberofcorrections             : Maximum number of corrections.
      --bound, --lbound, --ubound       : Parameter constraint settings.
  - Optionally, specify a grid position file using "--gridposition" to override the default grid setup.
  - Additional parameters mirror many of the affine registration options.

Example Command:
  3DRegBsplines \
    -f fixedImage.nii.gz \
    -m movingImage.nii.gz \
    -o registeredImage.nii.gz \
    -v deformationField.nii.gz \
    --gridresolution 50 \
    --numberofthreads 4 \
    -a 1.0 \
    -A 1.0 \
    -l 0 \
    -L 0 \
    -n 0 \
    -N 0 \
    -p 0.1 \
    -b 64 \
    --bsplinecaching true \
    --costfunctionconvergencefactor 1e12 \
    --projectedgradienttolerance 1e-5 \
    --numberofevaluations 500 \
    --numberofcorrections 5 \
    --bound 2 \
    --lbound 0 \
    --ubound 0

Again, run the executable with "--help" to view all available options and their descriptions.

------------------------------------------------------------
Notes:
- Replace file names (e.g., fixedImage.nii.gz, movingImage.nii.gz) with your actual image file names.
- Adjust the options based on your specific registration needs.

## Running with Docker
[erosmontin/mplus:latest](https://hub.docker.com/r/erosmontin/mplus)

Pull the Docker image:
```bash
docker pull erosmontin/mplus:latest
```
Run the Docker container:

```bash

docker run --rm -v $(pwd):/data erosmontin/mplus:latest 
```

## Running with Singularity on HPC

You can build and run the project using Singularity by utilizing either the provided `Singularity.def` file located in the root of this repository or directly from the Docker image.

### Option 1: Build from Singularity definition file

```bash
sudo singularity build mPlus.sif Singularity.def
```

### Option 2: Build directly from Docker image

You can build the Singularity container directly from the Docker Hub image:

```bash
sudo singularity build mPlus.sif docker://erosmontin/mplus:latest
```

### Running the Singularity container:

To run the executable directly:

```bash
singularity run mPlus.sif [arguments]
```

Or explicitly invoke the command:

```bash
singularity exec mPlus.sif 3DRegBsplines [arguments]
```

Replace `[arguments]` with appropriate options as detailed in the usage section above.


# Contributors
[*Dr. Eros Montin, PhD*](http://me.biodimensional.com)\
[![GitHub](https://img.shields.io/badge/GitHub-erosmontin-blue)](https://github.com/erosmontin)\
[![ORCID](https://img.shields.io/badge/ORCID-0000--0002--1773--0064-green)](https://orcid.org/0000-0002-1773-0064)\
[![Scopus](https://img.shields.io/badge/Scopus-35604121500-orange)](https://www.scopus.com/authid/detail.uri?authorId=35604121500)



## Cite Us
```
Montin E, Belfatto A, Bologna M, et al. A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Med Biol Eng Comput. 2020;58(4):843-855. doi:10.1007/s11517-019-02109-4
```
and

```
Cavatorta C, Meroni S, Montin E, et al. Retrospective study of late radiation-induced damages after focal radiotherapy for childhood brain tumors. PLoS One. 2021;16(2):e0247748. Published 2021 Feb 26. doi:10.1371/journal.pone.0247748
```


## Notes
1. if you build your ITK4 you need set to true ITKV3_COMPATIBILITY and the shared libraries
