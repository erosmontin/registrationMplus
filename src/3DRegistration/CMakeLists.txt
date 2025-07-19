cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(3DRegistration LANGUAGES CXX)

# Build each algorithm sub‐folder
add_subdirectory(3DRegBsplines/src)
add_subdirectory(3DRegAffine/src)

