cmake_minimum_required(VERSION 3.20)
project(Entangle)

set(CMAKE_CXX_STANDARD 20)

add_compile_options(-DMKL_ILP64 -m64 -Wl,--no-as-needed)


# Flags of sources, headers and libraries
include_directories(
        /home/shifeng/Codes/0.Codes/ITensor
        /home/shifeng/Codes/0.Codes/eigen3
        /home/shifeng/Codes/0.Codes/intel/oneapi/mkl/latest/include
)

link_directories(
        /home/shifeng/Codes/0.Codes/intel/oneapi/mkl/latest/lib/intel64
        /home/shifeng/Codes/0.Codes/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin
        /home/shifeng/Codes/0.Codes/ITensor/lib
)

set(MKL -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_rt -lmkl_core -liomp5 -lpthread)
set(ITENSOR -litensor)
set(SOURCE_FILES main.cpp)

add_executable(Entangle ${SOURCE_FILES})
target_link_libraries(Entangle ${MKL} ${ITENSOR})