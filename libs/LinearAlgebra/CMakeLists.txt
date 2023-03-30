find_package(PkgConfig REQUIRED)

pkg_search_module(PETSC PETSc>3.5.0 IMPORTED_TARGET REQUIRED)
message(STATUS "Found PETSc ${PETSC_VERSION}")

pkg_get_variable(PETSC_INCLUDE_DIR PETSC includedir)
pkg_get_variable(PETSC_LIB_DIR PETSC libdir)
message(STATUS "PETSc include directory: ${PETSC_INCLUDE_DIR}")
message(STATUS "PETSc library directory: ${PETSC_LIB_DIR}")

find_package(MPI REQUIRED)

# cmake-format: off
configure_library(NAME LinearAlgebra
                  SOURCE_FILES Vector.cpp Matrix.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES ""
                  INTERFACE_LINK_LIBRARIES "${PETSC_LIB_DIR}/libpetsc.dylib" ${PROJECT_NAME}::Utility)
# cmake-format: on

target_include_directories(${PROJECT_NAME}_LinearAlgebra SYSTEM PUBLIC ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

add_subdirectory(tests)
