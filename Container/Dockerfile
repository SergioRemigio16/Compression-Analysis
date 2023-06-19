# Use an official base image with a C++ compiler and Linux OS
FROM ubuntu:latest

# Set the working directory in the container
WORKDIR /app

# Install required packages and tools
RUN apt-get update && \
    apt-get install -y \
        build-essential \
        cmake \
        git \
        libboost-all-dev \
        vim \
        dos2unix \
        curl \
        libfftw3-dev

# Clone the repository and build zfp
RUN git clone https://github.com/LLNL/zfp.git /app/zfp && \
    cd /app/zfp && \
    mkdir build && \
    cd build && \
    cmake -DZFP_WITH_OPENMP=OFF .. && \
    cmake --build . --config Release && \
    ctest

# Create a new project directory and set it as the working directory
WORKDIR /app/my_project

# Copy the CMakeLists.txt and main.cpp files from the local machine to the Docker image
COPY src/CMakeLists.txt src/main.cpp src/Utilities.cpp src/Utilities.h src/CompressionDecompression.cpp src/CompressionDecompression.h ./

# Build the project
RUN mkdir build && \
    cd build && \
    cmake .. && \
    make

# Copy the build_script.sh file from the local machine to the Docker image and make it executable
COPY src/build_script.sh ./
RUN chmod +x build_script.sh
RUN dos2unix build_script.sh

# Set the entry point to the interactive shell
ENTRYPOINT ["/bin/bash"]
