cmake_minimum_required(VERSION 2.8.12)
project(SofaMiscForcefieldDev)

find_package(SofaAdvanced REQUIRED)


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


add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

set(SOFAMISCFORCEFIELDDEV_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${SOFAMISCFORCEFIELDDEV_INCLUDE_DIRS}>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include/${PROJECT_NAME}>")

target_link_libraries(${PROJECT_NAME} SofaBoundaryCondition SofaSimpleFem SofaBaseLinearSolver SofaUserInteraction)

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${HEADER_FILES}")
set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_MISC_FORCEFIELD_DEV")



## Install rules for the library and headers; CMake package configurations files
sofa_create_package(${PROJECT_NAME} 0.1 ${PROJECT_NAME} ${PROJECT_NAME})

install(DIRECTORY sofa
        DESTINATION include/${PROJECT_NAME}
        FILES_MATCHING PATTERN "*.h" PATTERN "*hpp" PATTERN "*.inl"
)
