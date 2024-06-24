# Registration based on MA+NGF+MSE


## Cite Us
```
Montin E, Belfatto A, Bologna M, et al. A multi-metric registration strategy for the alignment of longitudinal brain images in pediatric oncology. Med Biol Eng Comput. 2020;58(4):843-855. doi:10.1007/s11517-019-02109-4
```
and

```
Cavatorta C, Meroni S, Montin E, et al. Retrospective study of late radiation-induced damages after focal radiotherapy for childhood brain tumors. PLoS One. 2021;16(2):e0247748. Published 2021 Feb 26. doi:10.1371/journal.pone.0247748

```


## Installation

1. install the dependecnies 
    ```
    apt-get update && apt-get install -y cmake build-essential libinsighttoolkit4-dev
    ```
1. build the software in the root directory
    ```
    mkdir build && cd build
    cmake .. && make -j4 && make install
    ```

## Usage
1. Non rigid 3D registration
    ```
    bin/3DRecBsplines_new  -f FixedImage -m MovingImage -o outputRegisteredImage -v deformationfield  --numberofthreads 20 -p 01 -b 64 -B true --explicitPDFderivatives true -l 0.1 -L 0 -g 30 -N 1
    ```


## Docker version v1
https://hub.docker.com/r/erosmontin/mplus
erosmontin/mplus:tagname

## Docker version v2 coming