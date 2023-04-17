# cmake-format: off
configure_library(NAME Mesh
                  SOURCE_FILES Mesh.cpp Element.cpp Triangle.cpp Line.cpp Tetrahedron.cpp TriangleOrder2.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES Eigen3::Eigen
                  INTERFACE_LINK_LIBRARIES ${PROJECT_NAME}::Utility)
# cmake-format: on

add_subdirectory(tests)
