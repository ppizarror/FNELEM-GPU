# Make properties
cmake_minimum_required(VERSION 3.10)
project(FNELEM-GPU)
set(CMAKE_CXX_STANDARD 14)

# Math library
set(FNELEM_MATH
        fnelem/math/fematrix.cpp
        fnelem/math/fematrix_utils.cpp
        fnelem/math/matrix_inversion_cpu.cpp
        # fnelem/math/matrix_inversion_cuda.cu
        )

# Model library
set(FNELEM_MODEL
        fnelem/model/base/model_component.cpp
        fnelem/model/node/node.cpp
        )

# Define tests
add_executable(TEST-FEMATRIX test/math/test_fematrix.cpp ${FNELEM_MATH})
add_executable(TEST-FEMATRIX-UTILS test/math/test_fematrix_utils.cpp ${FNELEM_MATH})
add_executable(TEST-MODEL-COMPONENT test/model/base/test_model_component.cpp ${FNELEM_MODEL})
add_executable(TEST-NODE test/model/node/test_node.cpp ${FNELEM_MODEL})