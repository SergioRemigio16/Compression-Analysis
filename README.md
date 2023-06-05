# Compression-Analysis
Research of compression techniques for small matrices (3x7x7).

# Overview Slides
https://docs.google.com/presentation/d/1Pf6zKoRNxAkovl1baRmhLx5bBF0J83x5eqgKoctOmEM/edit?usp=drive_link

## Problem
Sending ghost points between CPUs using MPI is time-consuming.

## Solution 
Spend CPU cycles compressing/decompressing ghost points to decrease 
time spent transferring data. 

## Compression Techniques

ZFP library - Fast generic library designed for numerical 3D arrays that exploits patterns in floating-point data
SVD - Good compression if data has strong linear patterns
FFT - Transforms data to frequency domain if data has strong periodic components
Chebyshev - Good compression for smooth continuous functions

## Data
Physics simulation of black hole collisions in 3 dimensions 
using doubles and ran in linux using C++. The data is expected 
to emerge from wave functions. The data closer to the black 
holes is at a much higher resolution than to the rest of the 
simulated area to decrease computation time. This means points 
closer to the black holes are much closer together while points 
further from the black holes are more spaced out. There is 
much less activity away from the black holes so data is expected 
to be flat but still come from wave functions.

# Dockerfile notes
The DockerFile downloads any necessary packages and libraries 
to build and link the zfp library to the main.cpp program. 
The vimrc commands in the dockerfile are not necessary to build and run this 
program and can be removed if needed.

# How to build Docker container and run interactive mode
1. Open cmd at directory containing Dockerfile 
2. run the command "docker build -t zfp_testing ."
3. run the command "docker run -it zfp_testing"

# How to run program 
Once you run the docker container in interactive mode, 
you can use the following command to build and run the main.cpp program
"./build_script.sh"


# WindowsVisualStudioSolution
The WindowsVisualStudioSolution folder contains the same code
inside the Container except that it is used for windows development.
To successfully build and run this program in windows you need to install, 
build, and link the zfp library on your machine.
https://zfp.readthedocs.io/en/release1.0.0/installation.html