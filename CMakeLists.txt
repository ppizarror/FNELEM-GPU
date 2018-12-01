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

# ELEMENT LIBRARY
set(FNELEM_ELEMENT
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        ${FNELEM_NODE}
        fnelem/model/elements/element.cpp
        fnelem/model/elements/membrane.cpp
        )

# RESTRAINTS LIBRARY
set(FNELEM_RESTRAINTS
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        ${FNELEM_NODE}
        fnelem/model/restraints/node_restraint.cpp
        )

# LOADS LIBRARY
set(FNELEM_LOAD
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        ${FNELEM_NODE}
        fnelem/model/load/load.cpp
        )

# ALL LIBRARIES
set(FNELEM_ALL
        ${FNELEM_MATH}
        ${FNELEM_MODELCOMPONENT}
        ${FNELEM_NODE}
        ${FNELEM_ELEMENT}
        ${FNELEM_RESTRAINTS}
        )

# DEFINE TESTS
add_executable(TEST-ALL test/test_suite.cpp ${FNELEM_ALL})
add_executable(TEST-CUDA test/test_cuda.cpp ${FNELEM_CUDA})
add_executable(TEST-ELEMENTS test/model/elements/__elements__.cpp ${FNELEM_ELEMENT})
add_executable(TEST-FEMATRIX test/math/__math__.cpp ${FNELEM_MATH})
add_executable(TEST-LOAD test/model/load/__load__.cpp ${FNELEM_LOAD})
add_executable(TEST-MODEL-COMPONENT test/model/base/__base__.cpp ${FNELEM_MODELCOMPONENT})
add_executable(TEST-NODE test/model/node/__node__.cpp ${FNELEM_NODE})
add_executable(TEST-NODE-RESTRAINT test/model/restraints/__restraints__.cpp ${FNELEM_RESTRAINTS})