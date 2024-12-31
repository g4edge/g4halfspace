#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.21)
project(geant4halfspace)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available Vis drivers by default
# You can set WITH_GEANT4_VIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_VIS "Build example with Geant4 Vis drivers" ON)
if(WITH_GEANT4_VIS)
  find_package(Geant4 REQUIRED gdml vis_all)
else()
  find_package(Geant4 REQUIRED gdml)
endif()

#----------------------------------------------------------------------------
# Add external 3rd party libraries
#

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
find_package(CGAL 5 REQUIRED)
find_package(Catch2 REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})

file(GLOB HalfSpaceSources ${PROJECT_SOURCE_DIR}/src/*cc)
file(GLOB HalfSpaceHeaders ${PROJECT_SOURCE_DIR}/include/*.hh)
add_library(g4halfspace SHARED ${HalfSpaceSources} ${HalfSpcaeHeaders})
target_compile_options(g4halfspace PRIVATE -Werror )
target_link_libraries(g4halfspace ${Geant4_LIBRARIES} CGAL::CGAL Eigen3::Eigen)

file(GLOB HalfSpaceTestSources ${PROJECT_SOURCE_DIR}/test/*cc)
file(GLOB HalfSpaceTestHeaders ${PROJECT_SOURCE_DIR}/test/*hh)
add_executable(g4halfspace_test ${HalfSpaceTestSources} ${HalfSpaceTestHeaders})
target_compile_options(g4halfspace_test PRIVATE -Werror)
target_link_libraries(g4halfspace_test ${Geant4_LIBRARIES} g4halfspace)

# Copy over test files
set(TEST-INPUTS 00_trans.txt 01_aarbox.txt 02_aarbox_od.txt 03_rbox.txt 04_rbox_od.txt
                05_sphere.txt 06_cone.txt 07_cone_od.txt 08_ellipsoid.txt
                09_ellipsoid_od.txt 10_wedge.txt 11_wedge_od.txt 12_arbitrary.txt
                13_arbitrary_od.txt)
foreach(_input ${TEST-INPUTS})
  configure_file(
          ${PROJECT_SOURCE_DIR}/test-input/${_input}
          ${PROJECT_BINARY_DIR}/${_input}
          COPYONLY
  )
endforeach()

# Copy over macros
set(SCRIPTS vis.mac)
foreach(_script ${SCRIPTS})
  configure_file(
          ${PROJECT_SOURCE_DIR}/${_script}
          ${PROJECT_BINARY_DIR}/${_script}
          COPYONLY
  )
endforeach()

