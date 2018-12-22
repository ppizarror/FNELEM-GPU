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

# MODEL BASE LIBRARY
set(FNELEM_MODEL_BASE
        fnelem/model/base/model.cpp
        fnelem/model/base/model_component.cpp
        )

# MODEL NODE LIBRARY
set(FNELEM_MODEL_NODES
        ${FNELEM_MATH}
        ${FNELEM_MODEL_BASE}
        fnelem/model/nodes/node.cpp
        )

# MODEL ELEMENT LIBRARY
set(FNELEM_MODEL_ELEMENTS
        ${FNELEM_MATH}
        ${FNELEM_MODEL_BASE}
        ${FNELEM_MODEL_NODES}
        fnelem/model/elements/element.cpp
        fnelem/model/elements/membrane.cpp
        )

# MODEL RESTRAINTS LIBRARY
set(FNELEM_MODEL_RESTRAINTS
        ${FNELEM_MATH}
        ${FNELEM_MODEL_BASE}
        ${FNELEM_MODEL_NODES}
        fnelem/model/restraints/restraint.cpp
        fnelem/model/restraints/restraint_node.cpp
        )

# MODEL LOADS LIBRARY
set(FNELEM_MODEL_LOADS
        ${FNELEM_MATH}
        ${FNELEM_MODEL_BASE}
        ${FNELEM_MODEL_NODES}
        ${FNELEM_MODEL_ELEMENTS}
        fnelem/model/loads/load.cpp
        fnelem/model/loads/load_membrane_distributed.cpp
        fnelem/model/loads/load_node.cpp
        fnelem/model/loads/load_pattern.cpp
        fnelem/model/loads/load_pattern_constant.cpp
        )

# ALL LIBRARIES
set(FNELEM
        ${FNELEM_MATH}
        ${FNELEM_MODEL_BASE}
        ${FNELEM_MODEL_NODES}
        ${FNELEM_MODEL_ELEMENTS}
        ${FNELEM_MODEL_RESTRAINTS}
        ${FNELEM_MODEL_LOADS}
        )

# DEFINE TESTS
add_executable(TEST-ALL test/test_suite.cpp ${FNELEM})
add_executable(TEST-BASE test/model/base/__base__.cpp ${FNELEM})
add_executable(TEST-CUDA test/test_cuda.cpp ${FNELEM_CUDA})
add_executable(TEST-ELEMENTS test/model/elements/__elements__.cpp ${FNELEM_MODEL_ELEMENTS})
add_executable(TEST-FEMATRIX test/math/__math__.cpp ${FNELEM_MATH})
add_executable(TEST-LOADS test/model/loads/__loads__.cpp ${FNELEM_MODEL_LOADS})
add_executable(TEST-NODES test/model/nodes/__nodes__.cpp ${FNELEM_MODEL_NODES})
add_executable(TEST-RESTRAINTS test/model/restraints/__restraints__.cpp ${FNELEM_MODEL_RESTRAINTS})