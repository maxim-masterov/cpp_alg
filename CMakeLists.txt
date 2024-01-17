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
    set(CMAKE_CXX_FLAGS     "-O3 -g -march=core-avx2")
    if(NOT USE_STL)
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(OpenMP_CXX_FLAGS "-fopenmp")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
    # set(CMAKE_CXX_FLAGS  "-O3 -march=core-avx2 -stdpar=multicore --experimental-stdpar --gcc-toolchain=$ENV{EBROOTGCCCORE}/bin")
    set(CMAKE_CXX_FLAGS  "-O3 -g -stdpar=multicore --experimental-stdpar --gcc-toolchain=$ENV{EBROOTGCCCORE}/bin")
    set(OpenMP_CXX_FLAGS "")
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES Intel)
    set(CMAKE_CXX_FLAGS "-O3 -g -march=core-avx2 -qopenmp")
    if(NOT USE_STL)
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp")
        set(OpenMP_CXX_FLAGS "-qopenmp")
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(CMAKE_CXX_FLAGS "-O3 -g3 -march=native")
    set(CMAKE_C_FLAGS "-O3 -g3 -march=native")
    if(NOT USE_STL)
        # set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
        set(OpenMP_CXX_FLAGS "-fopenmp")
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

option(USE_STL "Use STL algorithms instead of OpenMP in the code" OFF) # Disabled by default

if(USE_STL)
    ADD_DEFINITIONS(-DUSE_STL)
    if(NOT CMAKE_CXX_COMPILER_ID MATCHES NVHPC)
        target_link_libraries(${CMAKE_PROJECT_NAME} tbb)
    endif()
else()
    target_link_libraries(${CMAKE_PROJECT_NAME})
endif(USE_STL)

message("")
message("//===================================================")
message("// ${PROJECT_NAME} build configuration:")
message("// ${PROJECT_NAME} version : ${PROJECT_VERSION}")
message("//===================================================")
message("  CMake version          : ${CMAKE_VERSION}")
message("  C++ Compiler           : ${CMAKE_CXX_COMPILER_ID}")
message("                           ${CMAKE_CXX_COMPILER_VERSION}")
message("                           ${CMAKE_CXX_COMPILER}")
message("  CXX STANDARD           : ${CMAKE_CXX_STANDARD}")
message("  CMAKE_CXX_FLAGS        : ${CMAKE_CXX_FLAGS}")
message("  OpenMP compile flags   : ${OpenMP_CXX_FLAGS}")
message("")