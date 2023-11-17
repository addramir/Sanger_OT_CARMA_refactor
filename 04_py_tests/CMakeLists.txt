cmake_minimum_required(VERSION 3.12)
project(ind_normal_carma)

# Specify the path to the ArmadilloConfig.cmake file
list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/include/")

set(CMAKE_OSX_ARCHITECTURES "x86_64")

find_package(Armadillo REQUIRED CONFIG)
if(NOT Armadillo_FOUND)
    message(FATAL_ERROR "Armadillo library not found.")
endif()



set(PYBIND11_CPP_STANDARD -std=c++17)
find_package(pybind11 REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)

# Set the correct Python include directory
set(PYTHON_INCLUDE_DIR ${Python_INCLUDE_DIRS})

pybind11_add_module(ind_normal_carma bindings.cpp)
target_include_directories(ind_normal_carma PRIVATE ${ARMADILLO_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIR})
target_link_libraries(ind_normal_carma PRIVATE ${ARMADILLO_LIBRARIES})
