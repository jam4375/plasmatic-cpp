# cmake-format: off
configure_library(NAME ProblemTypes
                  SOURCE_FILES HeatEq2D.cpp HeatEq3D.cpp Mechanical.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES Eigen3::Eigen
                  INTERFACE_LINK_LIBRARIES ${PROJECT_NAME}::Mesh ${PROJECT_NAME}::LinearAlgebra)
# cmake-format: on

add_subdirectory(tests)
