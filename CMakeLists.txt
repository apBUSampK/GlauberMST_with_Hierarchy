#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(GlauberMST)
set(CMAKE_PLATFORM_INDEPENDENT_CODE ON)

#----------------------------------------------------------------------------
# Find ROOT (required package)
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
find_package(ROOT)
find_package(ROOT CONFIG REQUIRED)
find_package(ROOT COMPONENTS MathMore)
include(${ROOT_USE_FILE})
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#

add_subdirectory(TGlauber)
link_directories(${ROOT_LIBRARY_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#

#add_library(TGlauNucleon SHARED ${PROJECT_SOURCE_DIR}/src/TGlauNucleon.cc ${PROJECT_SOURCE_DIR}/include/TGlauNucleon.hh)
#add_library(TGlauNucleus SHARED ${PROJECT_SOURCE_DIR}/src/TGlauNucleus.cc ${PROJECT_SOURCE_DIR}/include/TGlauNucleus.hh)
#add_library(TGlauberMC SHARED ${PROJECT_SOURCE_DIR}/src/TGlauberMC.cc ${PROJECT_SOURCE_DIR}/include/TGlauberMC.hh)
add_executable(GlauberMST GlauberMST.cc ${sources} ${headers})
#target_link_libraries(GRATE TGlauNucleon TGlauNucleus TGlauberMC ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} TGlauberLib)
target_link_libraries(GlauberMST ${ROOT_LIBRARIES} TGlauberLib)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS GlauberMST DESTINATION bin)

