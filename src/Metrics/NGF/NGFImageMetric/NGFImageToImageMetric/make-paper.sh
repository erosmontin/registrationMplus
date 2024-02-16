#!/bin/bash 
mkdir -p build && \
pushd build && \
cmake ..  && \
make -j 2  && \
mkdir Doc && \
cd Doc  && \ 
ln -s ../../Doc/* . && \
make && \
popd && \
ln -s build/Doc/NGFImageMetric.pdf . 

