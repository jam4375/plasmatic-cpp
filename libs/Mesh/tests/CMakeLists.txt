# cmake-format: off
configure_test_executable(NAME MeshTest
                          SOURCE_FILES main.cpp LineOrder2.cpp TriangleOrder2.cpp TetrahedronOrder2.cpp
                          SOURCE_DIR "."
                          BUILD_LINK_LIBRARIES ${PROJECT_NAME}::Mesh)
# cmake-format: on

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/assets/ DESTINATION ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/assets/Mesh)