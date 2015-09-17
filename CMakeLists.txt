cmake_minimum_required(VERSION 2.8.12)
project(SofaMiscForcefieldDev)

set(HEADER_FILES
    sofa/component/initMiscForcefieldDev.h
    sofa/component/container/ShapeMatchingRotationFinder.h
    sofa/component/container/ShapeMatchingRotationFinder.inl
    sofa/component/forcefield/ShapeMatchingForceField.h
    sofa/component/forcefield/ShapeMatchingForceField.inl
    sofa/component/forcefield/MappedBeamToTetraForceField.h
    sofa/component/forcefield/MappedBeamToTetraForceField.inl
    sofa/component/forcefield/Mapped3DoFForceField.h
    sofa/component/forcefield/Mapped3DoFForceField.inl
	sofa/component/controller/EvolutionParameterController.h
	sofa/component/controller/EvolutionParameterController.inl
)

set(SOURCE_FILES
    main.cpp
	sofa/component/initMiscForcefieldDev.cpp
	sofa/component/container/ShapeMatchingRotationFinder.cpp
    sofa/component/forcefield/ShapeMatchingForceField.cpp
    sofa/component/forcefield/MappedBeamToTetraForceField.cpp
    sofa/component/forcefield/Mapped3DoFForceField.cpp
	sofa/component/controller/EvolutionParameterController.cpp
)

set(README_FILE SofaMiscForcefieldDev.txt)

find_package(SofaAdvanced)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
add_definitions("-DSOFA_BUILD_MISC_FORCEFIELD_DEV")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")
target_link_libraries(${PROJECT_NAME} SofaBoundaryCondition SofaSimpleFem SofaBaseLinearSolver SofaUserInteraction)

## Install rules for the library and headers; CMake package configurations files
include(SofaPackageConfigHelpers)
sofa_create_package(SofaMiscForcefieldDev 0.1 ${PROJECT_NAME} SofaMiscForcefieldDev)
