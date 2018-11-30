cmake_minimum_required(VERSION 3.1)
project(JMM_CPPLibs)

set(JMM_CPPLibs_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/JMM_CPPLibs) #Changed 21/9/18 by Mirebeau (include->JMM_CPPLibs)

# Versionning
# ===========

file(STRINGS "${JMM_CPPLibs_INCLUDE_DIR}/JMM_CPPLibs_config.hpp" JMM_CPPLibs_version_defines
     REGEX "#define JMM_CPPLibs_VERSION_(MAJOR|MINOR|PATCH)")
foreach(ver ${JMM_CPPLibs_version_defines})
    if(ver MATCHES "#define JMM_CPPLibs_VERSION_(MAJOR|MINOR|PATCH) +([^ ]+)$")
        set(JMM_CPPLibs_VERSION_${CMAKE_MATCH_1} "${CMAKE_MATCH_2}" CACHE INTERNAL "")
    endif()
endforeach()
set(${PROJECT_NAME}_VERSION
    ${JMM_CPPLibs_VERSION_MAJOR}.${JMM_CPPLibs_VERSION_MINOR}.${JMM_CPPLibs_VERSION_PATCH})
message(STATUS "JMM_CPPLibs v${${PROJECT_NAME}_VERSION}")


add_library(JMM_CPPLibs INTERFACE)
target_include_directories(JMM_CPPLibs INTERFACE 
                           $<BUILD_INTERFACE:${JMM_CPPLibs_INCLUDE_DIR}>
                           $<INSTALL_INTERFACE:include>)

# Installation
# ============

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS JMM_CPPLibs EXPORT ${PROJECT_NAME}-targets) #Changed 21/9/18 by Mirebeau (include->JMM_CPPLibs) 

install(DIRECTORY ${JMM_CPPLibs_INCLUDE_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})

# Makes the project importable from the build directory
export(EXPORT ${PROJECT_NAME}-targets
       FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake")

set(JMM_CPPLibs_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for JMM_CPPLibsConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${JMM_CPPLibs_CMAKECONFIG_INSTALL_DIR})

# JMM_CPPLibs is header-only and does not depend on the architecture.
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${JMM_CPPLibs_CMAKECONFIG_INSTALL_DIR})
install(EXPORT ${PROJECT_NAME}-targets
        FILE ${PROJECT_NAME}Targets.cmake
DESTINATION ${JMM_CPPLibs_CMAKECONFIG_INSTALL_DIR})