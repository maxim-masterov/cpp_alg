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

add_executable(cpp_alg ${SOURCES})

#add_executable(CPP_alg_tmp main.cpp)
