# Make properties
cmake_minimum_required(VERSION 3.10)
project(FNELEM-GPU)
set(CMAKE_CXX_STANDARD 14)

# Define tests
add_executable(TEST-FEMATRIX test/math/test_fematrix.cpp)
add_executable(TEST-MODELCOMPONENT test/model/base/test_model_component.cpp)
add_executable(TEST-NODE test/model/node/test_node.cpp)