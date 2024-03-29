cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(SciFi)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})


find_package(Geant4 REQUIRED)
find_package(ROOT REQUIRED)


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
list(REMOVE_ITEM sources ${PROJECT_SOURCE_DIR}/src/._*)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
list(REMOVE_ITEM headers ${PROJECT_SOURCE_DIR}/include/._*)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(scifi scifi.cc ${sources} ${headers})
target_link_libraries(scifi ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
