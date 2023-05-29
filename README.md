# Compression-Analysis
Research of compression techniques for small matrices

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
To successfully build this solution you need to install, build, and link,
the zfp library on your machine.
https://zfp.readthedocs.io/en/release1.0.0/installation.html