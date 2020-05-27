cmake_minimum_required(VERSION 3.0)

# Set the project name
project(OOPS VERSION 1.1)

# Specify the C++ Standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Require Python
find_package(PythonInterp 3.0 REQUIRED)


# Set a default build type.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to Debug as none was specified.")
  set(CMAKE_BUILD_TYPE "Debug" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

option(USE_SDF "Enable exporting to SDF file" OFF)
option(BUILD_TESTS "Compile OOPS unit tests" OFF)

# Configure a header file that will pass on our settings.
configure_file(${PROJECT_SOURCE_DIR}/include/oopsconfig.h.in ${PROJECT_SOURCE_DIR}/include/oopsconfig.h)

# A macro to do a recursive search for subdirectories with CMakeLists.txt files.
macro(SUBDIRLIST result curdir)
  file(GLOB children RELATIVE ${curdir} ${curdir}/*)
  set(dirlist "")
  foreach(child ${children})
    if(IS_DIRECTORY ${curdir}/${child} AND EXISTS ${curdir}/${child}/CMakeLists.txt)
      list(APPEND dirlist ${child})
    endif()
  endforeach()
  set(${result} ${dirlist})
endmacro()

set(EXTRA_LIBS
   )

if(USE_SDF)
  list(APPEND EXTRA_LIBS bbhutil)
endif()

set(OOPS_INCLUDE_FILES
    include/grid.h
    include/domain.h
    include/rk4.h
    include/solver.h
    include/output.h
    include/cubic.h
    include/odedata.h
    include/solverdata.h
    include/ode.h
    include/parameters.h
    include/interpolator.h
    include/cubicinterpolator.h
    include/polynomialinterpolator.h
    include/paramreader.h
    include/fieldmap.h
   )

set(OOPS_SOURCE_FILES
    src/grid.cpp
    src/domain.cpp
    src/rk4.cpp
    src/ode.cpp
    src/odedata.cpp
    src/solverdata.cpp
    src/cubic.cpp
    src/interpolator.cpp
    src/cubicinterpolator.cpp
    src/polynomialinterpolator.cpp
    src/paramreader.cpp
    src/output.cpp
    src/fieldmap.cpp
   )

add_library(oops ${OOPS_INCLUDE_FILES} ${OOPS_SOURCE_FILES})
target_include_directories(oops PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# set(FLUID_INCLUDE_FILES
#     fluid/include/reconstruction.h
#     fluid/include/minmod.h
#     fluid/include/norecon.h
#     fluid/include/fluid.h
#     fluid/include/root.h
#     fluid/include/fluidparameters.h
#     fluid/include/fluidparser.h
#    )
# set(TEST_SOURCE_FILES
#     fluid/src/reconstruction.cpp
#     fluid/src/minmod.cpp
#     fluid/src/fluidmain.cpp
#     fluid/src/norecon.cpp
#     fluid/src/fluid.cpp
#     fluid/src/root.cpp
#     fluid/src/fluidparser.cpp
#    )
# 
# set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES} ${FLUID_INCLUDE_FILES})
# add_executable(Fluid ${SOURCE_FILES})
# target_link_libraries(Fluid ${EXTRA_LIBS})
# target_include_directories(Fluid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
#                            ${CMAKE_CURRENT_SOURCE_DIR}/fluid/include ${PROJECT_SOURCE_DIR})

set(WAVE_INCLUDE_FILES
    simplewave/include/wave.h
    )
set(TEST_SOURCE_FILES
    simplewave/src/wave.cpp
    simplewave/src/main.cpp
    )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES} ${WAVE_INCLUDE_FILES})
add_executable(Wave ${SOURCE_FILES})
target_link_libraries(Wave ${EXTRA_LIBS})
target_include_directories(Wave PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                            ${CMAKE_CURRENT_SOURCE_DIR}/simplewave/include ${PROJECT_SOURCE_DIR})

# Automatically add subdirectories containing CMakeLists.txt files and compile them.
SUBDIRLIST(SUBDIRS ${CMAKE_CURRENT_SOURCE_DIR})
foreach(subdir ${SUBDIRS})
  add_subdirectory(${subdir})
endforeach()
