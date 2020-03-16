#!/bin/bash

# Creating build directory
mkdir build
cd build

# Setting up environment
module load openmpi

# Building h5ToRoot example
cmake ..
make
