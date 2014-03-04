cmake_minimum_required(VERSION 2.8)

project("SofaMiscForcefieldDev")

include(${SOFA_CMAKE_DIR}/preProject.cmake)

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

set(README_FILE "SofaMiscForcefieldDev.txt")

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})

AddCompilerDefinitions("SOFA_BUILD_MISC_FORCEFIELD_DEV")
AddLinkerDependencies(BoundaryCondition SimpleFem BaseLinearSolver UserInteraction )

include(${SOFA_CMAKE_DIR}/postProject.cmake)