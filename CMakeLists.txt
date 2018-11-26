# PROJECT PROPERTIES
project(FNELEM-GPU)

# MAKE PROPERTIES
cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)

# CUDA SOURCES
set(FNELEM_CUDA
        fnelem/math/matrix_inversion_cuda.cu
        main.cu
        )

# MATH LIBRARY
set(FNELEM_MATH
        fnelem/math/fematrix.cpp
        fnelem/math/fematrix_utils.cpp
        fnelem/math/matrix_inversion_cpu.cpp
        )

# MODEL LIBRARY
set(FNELEM_MODELCOMPONENT
        fnelem/model/base/model_component.cpp
        )

# NODE LIBRARY
set(FNELEM_NODE
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        fnelem/model/node/node.cpp
        )

# ALL LIBRARIES
set(FNELEM_ALL
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        ${FNELEM_NODE}
        )

# DEFINE TESTS
add_executable(TEST-ALL test/test_suite.cpp ${FNELEM_ALL})
add_executable(TEST-CUDA test/test_cuda.cpp ${FNELEM_CUDA})
add_executable(TEST-FEMATRIX test/math/test_fematrix.cpp ${FNELEM_MATH})
add_executable(TEST-FEMATRIX-UTILS test/math/test_fematrix_utils.cpp ${FNELEM_MATH})
add_executable(TEST-MODEL-COMPONENT test/model/base/test_model_component.cpp ${FNELEM_MODELCOMPONENT})
add_executable(TEST-NODE test/model/node/test_node.cpp ${FNELEM_NODE})