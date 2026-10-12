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

cmake_minimum_required(VERSION 3.15)
#option(USE_TRACY "Search for other modules in sibling directories?" OFF)
option(USE_EASYPROFILER "Search for other modules in sibling directories?" OFF)
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
set( USE_OPENMP_DEFAULT OFF ) # Use of OpenMP is considered experimental
set(HAVE_EASYPROFILER 0)
if(USE_EASYPROFILER)
	set(easy_profiler_DIR
	"${PROJECT_SOURCE_DIR}/thirdparty/easy_profiler/lib/cmake/easy_profiler/"
	CACHE LOCATION
	"Build and install easy_profiler to this directory")
	set(HAVE_EASYPROFILER 1)
	find_package(easy_profiler REQUIRED)
endif()
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

include(OpmInit)
OpmSetPolicies()

if(USE_MPI)
	set(HDF5_PREFER_PARALLEL TRUE)
endif()

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
# with the find module.
# Note: OpmLibMain (included below) pulls in "${project}-prereqs" at the
# correct point, so we must NOT include it manually here.
#include ("${project}-prereqs")


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
	if(MPI_FOUND AND HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
		message(WARNING "When building parallel OPM flow we need a "
			"parallel version of hdf5, but found only a serial one. "
			"Please install a parallel hdf5 library for MPI "
			"(e.g with apt-get install libhdf5-mpi-dev) and do a clean "
			"rebuild (build after \"make clean\"). Continuing with "
			"only normal restart without hdf5 file support.")
		set(HDF5_FOUND OFF)
		unset(HAVE_HDF5)
	endif()
	include (CMakeLists_files.cmake)
endmacro (files_hook)

macro (tests_hook)
endmacro (tests_hook)


# Deck input support (HAVE_ECL_INPUT) is unconditionally available in current opm-common;
# the (removed) HAVE_ECL_INPUT switch is still referenced by the local copy of
# MaterialLawManagerTable here, so define it for every target in this module.
add_compile_definitions(HAVE_ECL_INPUT=1)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

# OpmLibMain calls the modern, project-prefixed hook
# (opm-flowexperimental_prereqs_hook); the legacy unprefixed prereqs_hook above
# is never invoked, so wire the library up to its dependencies explicitly here.
# Linking opmsimulators transitively brings the opm-common / opm-grid / dune
# include directories and compile definitions needed by the sources.
target_link_libraries(opmflowexperimental PUBLIC opmsimulators)

include(OpmBashCompletion)

add_library(moduleVersionExp OBJECT opm/simulators/utils/moduleVersion.cpp)
set_property(TARGET moduleVersionExp PROPERTY POSITION_INDEPENDENT_CODE ON)
# Inherit include dirs / defines from opmsimulators so moduleVersion.cpp can
# find its header (this is a local copy of the opm-simulators source).
target_link_libraries(moduleVersionExp PRIVATE opmsimulators)

# Strictly we only depend on the update-version target,
# but this is not exposed in a super-build.
#add_dependencies(moduleVersion opmsimulators)
set(CMAKE_PROJECT_NAME "${PROJECT_NAME}")
add_custom_target(opm-flowexperimental_prepare)

opm_add_test(matrix_merger_example
		ONLY_COMPILE	
		ALWAYS_ENABLE	
		SOURCES
		examples/matrix_merger_example.cpp
		$<TARGET_OBJECTS:moduleVersionExp>
)
	
# Disabled: needs opm/simulators/linalg/system/ISTLSolverExperiment.hpp, whose local
# copy predates the system-solver move to opm-simulators and redefines the
# SystemTypes.hpp aliases (RRMatrix, SystemVector, ...) as non-templates.
#opm_add_test(flow_blackoil_tracy
#	ONLY_COMPILE
#	ALWAYS_ENABLE
#	DEPENDS opmsimulators
#	LIBRARIES opmsimulators
#	SOURCES
#	examples/flow_blackoil_tracy.cpp #opm/simulators/linalg/system/ISTLSolverExperiment.cpp opm/simulators/linalg/system/SystemPreconditioner.cpp
#	$<TARGET_OBJECTS:moduleVersionExp>
#	)
if(USE_EASYPROFILER)
	opm_add_test(flow_blackoil_easyprofiler
		ONLY_COMPILE
		ALWAYS_ENABLE
		DEPENDS opmsimulators
		LIBRARIES opmsimulators
		SOURCES
		examples/flow_blackoil_easyprofiler.cpp
		$<TARGET_OBJECTS:moduleVersionExp>
		)
endif()


opm_add_test(flow_blackoil_noelemctx
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_noelemctx.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_noelemctx PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(flow_blackoil_nocache
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_nocache.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_nocache PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)


opm_add_test(flow_blackoil_fast
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_fast.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_fast PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

# A/B control for flow_blackoil_fast: identical except IntensiveQuantities = master's.
opm_add_test(flow_blackoil_masterIQ
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_masterIQ.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_masterIQ PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(flow_blackoil_drygas
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_drygas.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

opm_add_test(flow_blackoil_fast_oilgas
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_fast_oilgas.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_fast_oilgas PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(flow_blackoil_fast_oilgas_energy
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_fast_oilgas_energy.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)
target_sources(flow_blackoil_fast_oilgas_energy PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)


opm_add_test(flow_blackoil_fast_wateroil
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_fast_wateroil.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)	

target_sources(flow_blackoil_fast_wateroil
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_propertyevaluation
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_propertyevaluation.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_propertyevaluation
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_eclpropertyevaluation
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_eclpropertyevaluation.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_eclpropertyevaluation
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_eclpropertyevaluation_lowlevel
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_eclpropertyevaluation_lowlevel.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_eclpropertyevaluation_lowlevel
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_relpermevaluation
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_relpermevaluation.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_relpermevaluation
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)



# Control for flow_blackoil_baked: same standard material-law manager, but
# generic relperm — isolates convergence effects of the manager swap vs baking.
opm_add_test(flow_blackoil_fast_stdlaw
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_fast_stdlaw.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

# Control for flow_blackoil_baked: master's Model (FIBlackOilModel) instead of
# BlackOilModelFvFast, isolating what the experimental model still contributes.
opm_add_test(flow_blackoil_baked_stdmodel
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_baked_stdmodel.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

# flow_blackoil_fast with baked relperm tables (FlowProblemSimpleBaked).
opm_add_test(flow_blackoil_baked
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/flow_blackoil_baked.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

# Standalone intensive-quantities benchmark: master IQ vs fast IQ vs
# region-batched passes on fabricated states from a real deck.
opm_add_test(bench_iq_evaluation
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/bench_iq_evaluation.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

opm_add_test(test_eclrelpermevaluation
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_eclrelpermevaluation.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_eclrelpermevaluation
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_eclrelpermevaluation_lowlevel
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_eclrelpermevaluation_lowlevel.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_eclrelpermevaluation_lowlevel
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)

opm_add_test(test_eclrelpermevaluation_flowdiag
	ONLY_COMPILE
	ALWAYS_ENABLE
	DEPENDS opmsimulators
	LIBRARIES opmsimulators
	SOURCES
	examples/material/test_eclrelpermevaluation_flowdiag.cpp
	$<TARGET_OBJECTS:moduleVersionExp>
	)

target_sources(test_eclrelpermevaluation_flowdiag
	PRIVATE opm/material/fluidmatrixinteractions/MaterialLawManagerTable.cpp)




set(HAVE_TRACY 0)
if(USE_TRACY)
	set(HAVE_TRACY 1)
	# Tracy is already linked (keyword signature) to every opm_add_test target
	# by opm_add_target_options -> use_tracy when USE_TRACY_PROFILER is ON, so
	# these explicit plain-signature links are redundant and would clash
	# (plain vs keyword target_link_libraries on the same target).
endif()

if(USE_EASYPROFILER)
	target_link_libraries(flow_blackoil_fast easyprofiler)
endif()

# opm_add_test(ebos
#		ONLY_COMPILE
#		ALWAYS_ENABLE
#		DEPENDS opmsimulators
#		LIBRARIES opmsimulators
#		SOURCES
#		examples/ebos.cpp
#		$<TARGET_OBJECTS:moduleVersionExp>
#		)


# opm_add_test(ebos_experiment
#		ONLY_COMPILE
#		ALWAYS_ENABLE
#		DEPENDS opmsimulators
#		LIBRARIES opmsimulators
#		SOURCES
#		examples/ebos_experiment.cpp
#		$<TARGET_OBJECTS:moduleVersionExp>
#		)


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
