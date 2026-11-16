cmake_minimum_required(VERSION 3.20)

sknano_add_analysis_module(
    NAME TauPolarization
    HEADERS
        include/TauPolarization/TauPolarimeter.h
        include/TauPolarization/GenTauBuilder.h
        include/TauPolarization/TauPolarization.h
    SOURCES
        src/TauPolarimeter.cc
        src/GenTauBuilder.cc
        src/TauPolarization.cc
    LINKDEF include/TauPolarization/LinkDef.hpp
    INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include
    ANALYZERS TauPolarization)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
