#File Generated by Hans J. Johnson
#Please contact hans-johnson@uiowa.edu for making enhancments/corrections
#PROJECT(NIFTI)
# Configure Dart testing support.
#INCLUDE(${CMAKE_ROOT}/Modules/Dart.cmake)

SUBDIRS(znzlib niftilib)

# add option to build shared libraries. The default is OFF to maintain the
# current build behavior
#OPTION ( BUILD_SHARED_LIBS "Toggle building shared libraries")

# Michael Hanke 2007-04-25 Note:
# 'Utilities' is merged into the include install path to be able to easily
# install into any other path if necessary. For the same reason lib, bin and
# include dir for standalone builds are only set if not yet defined.


INCLUDE_DIRECTORIES(niftilib)
INCLUDE_DIRECTORIES(znzlib)

SET(NIFTI_INSTALL_BIN_DIR /bin)
SET(NIFTI_INSTALL_INCLUDE_DIR /include)
SET(NIFTI_INSTALL_LIB_DIR /lib)

# #################### for zlib #############################
option(BuildOwnZlib "Build own zlib" on)
if ( BuildOwnZlib )

ADD_DEFINITIONS(-DHAVE_ZLIB)
SUBDIRS(zlib)
INCLUDE_DIRECTORIES(zlib)
INCLUDE_DIRECTORIES(${ZXHZLIB_SOURCE_DIR})
INCLUDE_DIRECTORIES(${ZXHZLIB_BINARY_DIR})

SET(NIFTI_ZLIB_LIBRARIES zlib ${NIFTI_ZLIB_LIBRARIES})

else ( BuildOwnZlib )
	# Find zlib
  FIND_PACKAGE(ZLIB)
  IF(ZLIB_FOUND)
  	ADD_DEFINITIONS(-DHAVE_ZLIB)
  	INCLUDE_DIRECTORIES(${ZLIB_INCLUDE_DIRS})
  	LINK_DIRECTORIES(${ZLIB_LIBRARIES})
  	SET(NIFTI_ZLIB_LIBRARIES ${ZLIB_LIBRARIES})
  else (ZLIB_FOUND) 
	MESSAGE("Didn't find zlib, build own zlib")
	SUBDIRS(zlib)
	INCLUDE_DIRECTORIES(zlib)
	set( BuildOwnZlib on )
  ENDIF(ZLIB_FOUND)
endif ( BuildOwnZlib )
# #################### for zlib #############################

IF(NOT DEFINED NIFTI_INSTALL_LIB_DIR)
  MESSAGE(FATAL_ERROR "Must set NIFTI_INSTALL_LIB_DIR.")
ENDIF(NOT DEFINED NIFTI_INSTALL_LIB_DIR)
IF(NOT DEFINED NIFTI_INSTALL_BIN_DIR)
  MESSAGE(FATAL_ERROR "Must set NIFTI_INSTALL_BIN_DIR.")
ENDIF(NOT DEFINED NIFTI_INSTALL_BIN_DIR)
IF(NOT DEFINED NIFTI_INSTALL_INCLUDE_DIR)
  MESSAGE(FATAL_ERROR "Must set NIFTI_INSTALL_INCLUDE_DIR.")
ENDIF(NOT DEFINED NIFTI_INSTALL_INCLUDE_DIR)







