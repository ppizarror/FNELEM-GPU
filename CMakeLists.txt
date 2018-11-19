# Make properties
cmake_minimum_required(VERSION 3.10)
project(FNELEM-GPU)
set(CMAKE_CXX_STANDARD 14)

# Define tests
add_executable(TEST-MODELCOMPONENT test/test_ModelComponent.cpp)
add_executable(TEST-NODE test/test_Node.cpp)