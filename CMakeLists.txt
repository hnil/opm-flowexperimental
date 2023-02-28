# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

# Mandatory call to project
project(opm-flowexperimental C CXX)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental

# The following was copied from CMakeLists.txt in opm-common.
# TODO: factor out the common parts in opm-common and opm-simulator as a cmake module

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)
find_package(opm-simulators REQUIRED)

include(OpmInit)
OpmSetPolicies()

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# Look for the opm-tests repository; if found the variable
# HAVE_OPM_TESTS will be set to true.
include(Findopm-tests)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include ("${project}-prereqs")


#find_package(amgcl)



# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
  opm_need_version_of ("dune-common")
  opm_need_version_of ("dune-istl")
  if(dune-fem_FOUND)
      opm_need_version_of ("dune-fem")
  endif()
  opm_need_version_of ("opm-models")
	opm_need_version_of ("opm-simulators")
  if(NOT fmt_FOUND)
    add_definitions(-DFMT_HEADER_ONLY)
    list(APPEND EXTRA_INCLUDES SYSTEM ${PROJECT_SOURCE_DIR}/external/fmtlib/include)
  endif()
  include_directories(${EXTRA_INCLUDES})

  include(UseDamaris)
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (fortran_hook)
endmacro (fortran_hook)

macro (files_hook)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)
 
    
# all setup common to the OPM library modules is done here
include (OpmLibMain)




include(OpmBashCompletion)


#add_library(moduleVersion OBJECT opm/simulators/utils/moduleVersion.cpp)
#set_property(TARGET moduleVersion PROPERTY POSITION_INDEPENDENT_CODE ON)

# Strictly we only depend on the update-version target,
# but this is not exposed in a super-build.
#add_dependencies(moduleVersion opmsimulators)
set(CMAKE_PROJECT_NAME "${PROJECT_NAME}")
add_custom_target(opm-flowexperimental_prepare)

opm_add_test(flow_blackoil_tracy
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  examples/flow_blackoil_tracy.cpp
  $<TARGET_OBJECTS:moduleVersion>
)

opm_add_test(flow_blackoil_noelemctx
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  examples/flow_blackoil_tracy.cpp
  $<TARGET_OBJECTS:moduleVersion>
)

opm_add_test(flow_blackoil_fast
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  examples/flow_blackoil_tracy.cpp
  $<TARGET_OBJECTS:moduleVersion>
)


opm_add_test(ebos
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  examples/ebos.cpp
  $<TARGET_OBJECTS:moduleVersion>
  )


opm_add_test(ebos_experiment
  ONLY_COMPILE
  ALWAYS_ENABLE
  DEPENDS opmsimulators
  LIBRARIES opmsimulators
  SOURCES
  examples/ebos_experiment.cpp
  $<TARGET_OBJECTS:moduleVersion>
  )


# if(dune-alugrid_FOUND)
#   if (NOT BUILD_FLOW_ALU_GRID)
#     set(FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF "FALSE")
#   else()
#     set(FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF "TRUE")
#   endif()

#   opm_add_test(flow_alugrid
#     ONLY_COMPILE
#     ALWAYS_ENABLE
#     DEFAULT_ENABLE_IF ${FLOW_ALUGRID_ONLY_DEFAULT_ENABLE_IF}
#     DEPENDS opmsimulators
#     LIBRARIES opmsimulators
#     SOURCES
#     flow/flow_blackoil_alugrid.cpp
#     $<TARGET_OBJECTS:moduleVersion>)
#   target_compile_definitions(flow_alugrid PRIVATE USE_ALUGRID)
# endif()



# must link libraries after target 'opmsimulators' has been defined
#install(DIRECTORY doc/man1 DESTINATION ${CMAKE_INSTALL_MANDIR}
#  FILES_MATCHING PATTERN "*.1")

