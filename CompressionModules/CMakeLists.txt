cmake_minimum_required(VERSION 3.10)
project(CompressionAnalysis)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)

# Include header files
include_directories(include)

# Include FFTW3 directories
include_directories(/usr/include)
link_directories(/usr/lib/x86_64-linux-gnu)

# Include ZFP directories
include_directories(/home/laptop/zfp/build)
link_directories(/home/laptop/zfp/build/lib)

# Include Blosc directories
include_directories(/home/sergio/Compression-Analysis/c-blosc/build/blosc_library/include)
link_directories(/home/sergio/Compression-Analysis/c-blosc/build/blosc_library/lib)

# Include Blosc directories
include_directories(${HOME}/Compression-Analysis/c-blosc/build/blosc_library/include)
link_directories(${HOME}/Compression-Analysis/c-blosc/build/blosc_library/lib)

# External libraries
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig)
pkg_check_modules(FFTW3 fftw3 REQUIRED)
find_package(ZFP REQUIRED)
find_package(OpenCV REQUIRED)  

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Generate the executable
add_executable(CompressionModules ${SOURCES})

# Link necessary libraries
target_link_libraries(CompressionModules
    PRIVATE
    Boost::boost
    Eigen3::Eigen
    ${FFTW3_LIBRARIES}
    ${ZFP_LIBRARIES}
    ${OpenCV_LIBS}  
    blosc  
)

# If FFTW3 include directories are needed
target_include_directories(CompressionModules PUBLIC ${FFTW3_INCLUDE_DIRS})

# If FFTW3 compile options are needed
target_compile_options(CompressionModules PUBLIC ${FFTW3_CFLAGS_OTHER})
