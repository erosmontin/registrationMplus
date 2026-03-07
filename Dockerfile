FROM ubuntu:jammy

# Install the necessary packages for ITK 5.x + C++17 build

RUN apt-get update && apt-get install -y --no-install-recommends \
    cmake build-essential \
    libinsighttoolkit5-dev \
    libpng-dev libjpeg-dev libtiff-dev libdcmtk-dev \
    libfltk1.3-dev libeigen3-dev \
    libboost-program-options-dev \
    python3 python3-pip \
    && rm -rf /var/lib/apt/lists/*

# Install Python package (pure-Python mode — no pybind11 compilation)
COPY python /python
RUN pip3 install --no-cache-dir /python

# Build C++ executables
COPY src /src
WORKDIR /bld
RUN cmake /src \
    -DCMAKE_BUILD_TYPE=Release \
    && make -j$(nproc)

ENV PATH="/bld/bin:${PATH}"

# Verify binaries exist after build
RUN ls -la /bld/bin/ || echo "Warning: /bld/bin not found; check cmake output above"

ENTRYPOINT ["/bld/bin/3DRegBsplines"]




