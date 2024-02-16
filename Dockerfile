from ubuntu:lunar

# Install the necessary packages for itk



RUN apt-get update && apt-get install -y cmake build-essential libinsighttoolkit4-dev


# build the software in the root directory
COPY src /src

WORKDIR /bld

RUN cmake /src && make -j4





