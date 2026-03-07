FROM ubuntu:jammy

# Install the necessary packages for itk



RUN apt-get update 
RUN apt-get install -y cmake build-essential libinsighttoolkit4-dev
RUN apt-get install libpng-dev libjpeg-dev libtiff-dev libdcmtk-dev libfltk1.3-dev libeigen3-dev -y
RUN apt-get install libboost-all-dev -y
COPY src /src

WORKDIR /bld

RUN cmake /src && make -j4

ENV PATH="/bld/bin:${PATH}"

# Verify binary exists after build
RUN ls -la /bld/bin/ || echo "Warning: /bld/bin not found; check cmake output above"

ENTRYPOINT ["/bld/bin/3DRegBsplines"]
#3DRegBsplines -m /tmp/dess.nii.gz -f /tmp/tse.nii -o /tmp/r.nii.gz -y 0 -g 15 -L 0 -l 0 -n 0 -N 0 -p 0.001 -v /tmp/VF.mhd --numberofthreads=8
#Option: NGFevaluator
#Value: 0
#




