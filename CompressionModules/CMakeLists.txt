cmake_minimum_required(VERSION 3.10)
project(CompressionAnalysis)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Include header files
include_directories(include)

# Include FFTW3 directories
include_directories(/usr/include)  # For FFTW3 header files (optional, if FFTW headers are not found)
link_directories(/usr/lib/x86_64-linux-gnu)  # For FFTW3 library files

# External libraries
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig)
pkg_check_modules(FFTW3 fftw3 REQUIRED)
find_package(ZFP REQUIRED)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Generate the executable
add_executable(CompressionModules ${SOURCES})

# Link necessary libraries
target_link_libraries(CompressionModules
    PRIVATE
    Boost::boost
    Eigen3::Eigen
    ${FFTW3_LIBRARIES}  # This is for FFTW3
    ${ZFP_LIBRARIES}
)

# If FFTW3 include directories are needed, you can include them this way:
target_include_directories(CompressionModules PRIVATE ${FFTW3_INCLUDE_DIRS})

# If FFTW3 compile options are needed, you can include them this way:
target_compile_options(CompressionModules PRIVATE ${FFTW3_CFLAGS_OTHER})