# Make properties
cmake_minimum_required(VERSION 3.10)
project(FNELEM-GPU)
set(CMAKE_CXX_STANDARD 14)

# Source files
set(FNELEM_MATH
        fnelem/math/fematrix.cpp
        fnelem/math/fematrix_utils.cpp
        fnelem/math/matrix_inversion_cpu.cpp
        # fnelem/math/matrix_inversion_cuda.cu
        )
set(FNELEM_MODEL_BASE
        )

# Define tests
add_executable(TEST-FEMATRIX test/math/test_fematrix.cpp ${FNELEM_MATH})
add_executable(TEST-FEMATRIX-UTILS test/math/test_fematrix_utils.cpp ${FNELEM_MATH})
add_executable(TEST-MODEL-COMPONENT test/model/base/test_model_component.cpp)
add_executable(TEST-NODE test/model/node/test_node.cpp)