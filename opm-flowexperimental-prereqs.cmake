# defines that must be present in config.h for our headers
set (opm-flowexperimental_CONFIG_VAR
  HAVE_OPM_GRID
  HAVE_PTHREAD
  HAVE_MPI
  HAVE_PETSC
  HAVE_SUITESPARSE_UMFPACK_H
  HAVE_DUNE_ISTL
  DUNE_ISTL_VERSION_MAJOR
  DUNE_ISTL_VERSION_MINOR
  DUNE_ISTL_VERSION_REVISION
  HAVE_SUITESPARSE_UMFPACK
  HAVE_HDF5
  HAVE_TRACY
  )

# dependencies, found via direct find_package calls (modern opm-common no
# longer ships the find_package_deps() helper / ${project}_DEPS mechanism).
find_package(Boost COMPONENTS date_time unit_test_framework REQUIRED)
find_package(dune-common REQUIRED)
find_package(dune-istl REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(MPI)
find_package(SuiteSparse COMPONENTS UMFPACK)
find_package(opm-common REQUIRED)
find_package(opm-grid REQUIRED)
find_package(opm-simulators REQUIRED)
find_package(HDF5)

# Note: have to look for alugrid before dune-fem.
find_package(dune-alugrid)
find_package(dune-fem)
