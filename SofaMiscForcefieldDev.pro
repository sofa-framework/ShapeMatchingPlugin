# Target is a plugin: SofaMiscForcefieldDev

load(sofa/pre)
defineAsPlugin(SofaMiscForcefieldDev)

TARGET = SofaMiscForcefieldDev

DEFINES += SOFA_BUILD_MISC_FORCEFIELD_DEV

README_FILE = SofaMiscForcefieldDev.txt

HEADERS += sofa/component/initMiscForcefieldDev.h \
           sofa/component/container/ShapeMatchingRotationFinder.h \
           sofa/component/container/ShapeMatchingRotationFinder.inl \
           sofa/component/forcefield/ShapeMatchingForceField.h \
           sofa/component/forcefield/ShapeMatchingForceField.inl \
           sofa/component/forcefield/MappedBeamToTetraForceField.h \
           sofa/component/forcefield/MappedBeamToTetraForceField.inl

SOURCES += main.cpp \
		   sofa/component/initMiscForcefieldDev.cpp \
		   sofa/component/container/ShapeMatchingRotationFinder.cpp \
           sofa/component/forcefield/ShapeMatchingForceField.cpp \
           sofa/component/forcefield/MappedBeamToTetraForceField.cpp


# Make sure there are no cross-dependencies
INCLUDEPATH -= $$SOFA_INSTALL_INC_DIR/applications
DEPENDPATH -= $$SOFA_INSTALL_INC_DIR/applications

#exists(component-local.cfg): include(component-local.cfg)

load(sofa/post)
 
