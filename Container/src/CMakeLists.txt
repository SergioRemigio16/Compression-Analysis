cmake_minimum_required(VERSION 3.15)

project(my_project VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Locate the ZFP and FFTW3 libraries.
list(APPEND CMAKE_PREFIX_PATH "/app/zfp/build")
find_package(ZFP REQUIRED)
find_library(FFTW3_LIB fftw3)

add_executable(my_project_executable
    main.cpp
    Utilities.cpp
    CompressionDecompression.cpp
)

# Include ZFP headers and link ZFP and FFTW3 libraries to your executable.
target_include_directories(my_project_executable PRIVATE ${ZFP_INCLUDE_DIRS})
target_link_libraries(my_project_executable PRIVATE ${ZFP_LIBRARIES} ${FFTW3_LIB})
