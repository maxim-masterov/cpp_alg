cmake_minimum_required(VERSION 3.23)
project(cpp_alg)

set(CMAKE_CXX_STANDARD 23)

file(GLOB SOURCES
        DataTypes/field.h
        DataTypes/matrix.h
        DataTypes/vector.h
        General/dimensions.h
        General/helpers.h
        General/macro.h
        General/structs.h
        General/termination.h
        IO/io.h
        MPI/common.h
        MPI/Decomposition/decomposition.h
        Solver/solver.h
        System/system.h
        DataTypes/field.cpp
        DataTypes/matrix.cpp
        DataTypes/vector.cpp
        General/dimensions.cpp
        General/helpers.cpp
        IO/io.cpp
        MPI/common.cpp
        MPI/Decomposition/decomposition.cpp
        Solver/solver.cpp
        System/system.cpp
        main.cpp
    )

add_executable(${CMAKE_PROJECT_NAME} ${SOURCES})

if(CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(CMAKE_CXX_FLAGS         "-O3 -march=core-avx2 -fopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=core-avx2 -fopenmp")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    set(CMAKE_CXX_FLAGS        "-O3 -march=core-avx2 -stdpar=multicore --experimental-stdpar --std=c++20 --gcc-toolchain=$ENV{EBROOTGCCCORE}/bin")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=core-avx2 -stdpar=multicore --experimental-stdpar --std=c++20 --gcc-toolchain=$ENV{EBROOTGCCCORE}/bin")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS         "-O3 -g -march=core-avx2 -qopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g -qopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=core-avx2 -qopenmp")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS         "-O3 -g3 -march=native -fopenmp")
    set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -fopenmp")
endif()

option(USE_STL "Use STL algorithms instead of OpenMP in the code" OFF) # Disabled by default

if(USE_STL)
    ADD_DEFINITIONS(-DUSE_STL)
    target_link_libraries(${CMAKE_PROJECT_NAME} tbb)
else()
    target_link_libraries(${CMAKE_PROJECT_NAME})
endif(USE_STL)

# add_definitions(-DUSE_STL)
# target_link_libraries(cpp_alg tbb)

# target_link_libraries(cpp_alg)
