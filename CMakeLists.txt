cmake_minimum_required(VERSION 3.12)
project(ShapeMatchingPlugin LANGUAGES CXX)

set(SHAPEMATCHINGPLUGIN_SRC "src/${PROJECT_NAME}")

find_package(Sofa.Core REQUIRED)

set(HEADER_FILES
    ${SHAPEMATCHINGPLUGIN_SRC}/config.h.in
    ${SHAPEMATCHINGPLUGIN_SRC}/initShapeMatchingPlugin.h
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingRotationFinder.h
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingRotationFinder.inl
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingForceField.h
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingForceField.inl
)

set(SOURCE_FILES
	${SHAPEMATCHINGPLUGIN_SRC}/initShapeMatchingPlugin.cpp
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingRotationFinder.cpp
    ${SHAPEMATCHINGPLUGIN_SRC}/ShapeMatchingForceField.cpp
)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC Sofa.Core)

sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${Sofa_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR "${PROJECT_NAME}"
    RELOCATABLE "plugins"
)
