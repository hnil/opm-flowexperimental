# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
# MAIN_SOURCE_FILES     List of compilation units which will be included in
#                       the library. If it isn't on this list, it won't be
#                       part of the library. Please try to keep it sorted to
#                       maintain sanity.
#
# TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
# TEST_DATA_FILES       Files from the source three that should be made
#                       available in the corresponding location in the build
#                       tree in order to run tests there.
#
# EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#                       build, but which is not part of the library nor is
#                       run as tests.
#
# PUBLIC_HEADER_FILES   List of public header files that should be
#                       distributed together with the library. The source
#                       files can of course include other files than these;
#                       you should only add to this list if the *user* of
#                       the library needs it.

# originally generated with the command:
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/material/fluidmatrixinteractions/EclMaterialLawManagerTable.cpp
	)


# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort



# originally generated with the command:
# find tutorials examples -name '*.c*' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
#     examples/flow_test.cpp
#     examples/flow_test2.cpp
#     examples/flow_test_2ph.cpp
	)
list (APPEND  PROGRAM_SOURCE_FILES
 #    examples/flow_test.cpp
 #    examples/flow_test2.cpp
 #    examples/flow_test_2ph.cpp
		)
# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
		#eclproblemsimple.hh
	)

list(APPEND TEST_SOURCE_FILES
		# tests/material/test_eclpropertyevaluation.cpp
		# tests/material/test_eclpropertyevaluation_lowlevel.cpp
		# tests/material/test_eclpropertyevaluation_lowlevel_spe9.cpp
		# tests/material/test_relpermevaluation.cpp
		# tests/material/test_eclrelpermevaluation.cpp
		# tests/material/test_eclrelpermevaluation_lowlevel.cpp
)
