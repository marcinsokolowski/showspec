cmake_minimum_required (VERSION 2.6)
project (showspec_standalone CXX Fortran)

# The version number.
set (showspec_standalone_VERSION_MAJOR 1)
set (showspec_standalone_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
# configure_file (
#  "${PROJECT_SOURCE_DIR}//antenna_pattern.h.in"
#  "${PROJECT_BINARY_DIR}/antenna_pattern.h"
#  )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")


# Set the output folder where your program will be created
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# include_directories(${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-O3 -D_UNIX -D_NO_ROOT_")

# HEALPIX:
# My local healpix :
SET(MYLOCAL_HEALPIX_INCLUDEPATH /opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/include/)
SET(MYLOCAL_HEALPIX_LIBPATH /opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/lib/)

FIND_PATH(HEALPIX_INCLUDE_DIR
NAMES
        healpix_base.h
PATHS
        /usr/include
        /usr/local/include
        ${MYLOCAL_HEALPIX_INCLUDEPATH}
)

FIND_LIBRARY(HEALPIX_LIBRARIES
NAMES
        healpix_cxx
PATHS
        /usr/lib
        /usr/local/lib
        ${MYLOCAL_HEALPIX_LIBPATH}
)
FIND_LIBRARY(HEALPIX_CXXSUPPORT_LIBRARIES
NAMES
        cxxsupport
PATHS
        /usr/lib
        /usr/local/lib
        ${MYLOCAL_HEALPIX_LIBPATH}
)

SET(HEALPIX_FOUND FALSE)

# IF(HEALPIX_INCLUDE_DIR AND HEALPIX_LIBRARIES AND HEALPIX_CXXSUPPORT_LIBRARIES)
IF(HEALPIX_INCLUDE_DIR)
   SET(HEALPIX_LIBRARIES ${HEALPIX_LIBRARIES}
   ${HEALPIX_CXXSUPPORT_LIBRARIES}
   ${HEALPIX_SHARP_LIBRARIES}
   ${HEALPIX_FFTPACK_LIBRARIES}
   ${HEALPIX_CUTILS_LIBRARIES})
   SET(HEALPIX_FOUND TRUE)

   MESSAGE(STATUS "HealPix Found!")
   MESSAGE(STATUS "HEALPIX_INCLUDE_DIR=${HEALPIX_INCLUDE_DIR}")
   MESSAGE(STATUS "HEALPIX_LIBRARIES=${HEALPIX_LIBRARIES}")

   include_directories(${HEALPIX_INCLUDE_DIR})
ENDIF()

# ROOT OLD:
# if(EXISTS $ENV{ROOTSYS}/etc/cmake/)
#   MESSAGE(STATUS "ROOTSYS/etc/cmake/ found")
#   set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake/)
#   include($ENV{ROOTSYS}/etc/cmake/FindROOT.cmake)
#   find_package(ROOT)
#   include_directories(${ROOT_INCLUDE_DIR})
#else(EXISTS $ENV{ROOTSYS}/etc/cmake/)
#   MESSAGE(STATUS "ROOTSYS/etc/cmake/ not found")
#   include_directories($ENV{ROOTSYS}/include)
# endif(EXISTS $ENV{ROOTSYS}/etc/cmake/)

# libnova, fitsio etc :
find_library(NOVA_LIB NAMES nova)
find_library(FITSIO_LIB cfitsio REQUIRED)

add_executable(showspec_standalone showmap.cpp antenna_pattern.cpp user_defined_antenna_pattern.cpp skymap_cache.cpp ce300e.cpp showspec_math.cpp MollweideSkyMap.cpp bg_components.cpp bg_array.cpp bg_fits.cpp bg_geo.cpp libnova_interface.cpp bg_total_power.cpp bg_units.cpp bg_globals.cpp  bg_date.cpp weather_station.cpp bg_bedlam.cpp cvalue_vector.cpp bg_stat.cpp myparser.cpp mystring.cpp basestring.cpp myfile.cpp mysafekeytab.cpp cexcp.cpp mycmnglobals.cpp mydate.cpp mystrtable.cpp mykeytab.cpp myfract.cpp myenv.cpp basestructs.cpp random.cpp mylock.cpp mydate2.cpp myfits.cpp mathfunc.cpp mwa_beam_interface.cpp )
# target_link_libraries(showspec_standalone ${NOVA_LIB} ${FITSIO_LIB} ${ROOT_LIBRARIES} ${HEALPIX_LIBRARIES})
target_link_libraries(showspec_standalone ${NOVA_LIB} ${FITSIO_LIB} ${HEALPIX_LIBRARIES})

find_package(HDF5 COMPONENTS CXX REQUIRED)
if(HDF5_FOUND)
   MESSAGE(STATUS "HDF5 package found -> building version for MWA beam")
   set(CMAKE_CXX_FLAGS "-O3 -D_UNIX -D_NO_ROOT_")
   include_directories(${HDF5_CXX_INCLUDE_DIRS} ${HDF5_CXX_INCLUDE_DIR})
   SET_SOURCE_FILES_PROPERTIES( system.cpp PROPERTIES COMPILE_FLAGS -std=gnu++11 )
   add_executable(showspec_standalone_mwa showmap.cpp antenna_pattern.cpp user_defined_antenna_pattern.cpp skymap_cache.cpp ce300e.cpp showspec_math.cpp MollweideSkyMap.cpp bg_components.cpp bg_array.cpp bg_fits.cpp bg_geo.cpp libnova_interface.cpp bg_total_power.cpp bg_units.cpp bg_globals.cpp  bg_date.cpp weather_station.cpp bg_bedlam.cpp cvalue_vector.cpp bg_stat.cpp myparser.cpp mystring.cpp basestring.cpp myfile.cpp mysafekeytab.cpp cexcp.cpp mycmnglobals.cpp mydate.cpp mystrtable.cpp mykeytab.cpp myfract.cpp myenv.cpp basestructs.cpp random.cpp mylock.cpp mydate2.cpp myfits.cpp mathfunc.cpp mwa_beam_interface.cpp beam2016implementation.cpp system.cpp)
   # target_link_libraries(showspec_standalone_mwa ${NOVA_LIB} ${FITSIO_LIB} ${ROOT_LIBRARIES} ${HEALPIX_LIBRARIES})
   set(HDF5_LIB hdf5_cpp hdf5_serial)
   find_package(Boost 1.48 COMPONENTS filesystem system thread REQUIRED)
   target_link_libraries(showspec_standalone_mwa ${NOVA_LIB} ${FITSIO_LIB} ${HEALPIX_LIBRARIES} ${HDF5_C_LIBRARY} ${HDF5_LIB} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY})
   set_target_properties(showspec_standalone_mwa PROPERTIES COMPILE_FLAGS "-D_MWA_2016_BEAM_MODEL_")
else(HDF5_FOUND)
   MESSAGE(ERROR "HDF5 package not found - cannot build program showspec_standalone_mwa to integrate mwa beam with the GSM model")
endif(HDF5_FOUND)

# GSM model generator :
set(FC gfortran)
set(CMAKE_Fortran_FLAGS "-ffixed-line-length-200")
add_executable( skybkg_angelica_range gsm_range.f)

# SHOW ALL VARIABLES :
   get_cmake_property(_variableNames VARIABLES)
   foreach (_variableName ${_variableNames})
      message(STATUS "${_variableName}=${${_variableName}}")
   endforeach()


