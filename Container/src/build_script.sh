#!/bin/bash

# Change to the project directory
cd /app/my_project

# Create the build directory if it doesn't exist
mkdir -p build

# Change to the build directory
cd build

# Run CMake to configure the project
cmake ..

# Build the executable
make

# Run the executable
./my_project_executable

# Return to the project directory
cd ..

# Exit the script
exit 0
