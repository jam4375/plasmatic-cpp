find_package(PkgConfig REQUIRED)

pkg_search_module(PETSc PETSc>3.5.0 IMPORTED_TARGET REQUIRED)
if(${PETSc_VERSION})
    message(STATUS "Found PETSc ${PETSc_VERSION}")
endif()

pkg_get_variable(PETSc_INCLUDE_DIR PETSc includedir)
pkg_get_variable(PETSc_LIB_DIR PETSc libdir)
message(STATUS "PETSc include directory: ${PETSc_INCLUDE_DIR}")
message(STATUS "PETSc library directory: ${PETSc_LIB_DIR}")

find_package(MPI REQUIRED)

if(UNIX AND NOT APPLE)
    set(PETSC_LIB "${PETSc_LIB_DIR}/libpetsc.so")
endif()

if(APPLE)
    set(PETSC_LIB "${PETSc_LIB_DIR}/libpetsc.dylib")
endif()

# cmake-format: off
configure_library(NAME LinearAlgebra
                  SOURCE_FILES Vector.cpp Matrix.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES ""
                  INTERFACE_LINK_LIBRARIES ${PETSC_LIB} ${PROJECT_NAME}::Utility)
# cmake-format: on

target_include_directories(${PROJECT_NAME}_LinearAlgebra SYSTEM PUBLIC ${PETSC_INCLUDE_DIR} ${MPI_INCLUDE_PATH})

add_subdirectory(tests)
