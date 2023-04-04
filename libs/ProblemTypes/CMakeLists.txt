# cmake-format: off
configure_library(NAME ProblemTypes
                  SOURCE_FILES ProblemTypes.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES
                  INTERFACE_LINK_LIBRARIES ${PROJECT_NAME}::Mesh ${PROJECT_NAME}::LinearAlgebra)
# cmake-format: on

add_subdirectory(tests)
