# cmake-format: off
configure_test_executable(NAME MeshTest
                          SOURCE_FILES main.cpp
                          SOURCE_DIR "."
                          BUILD_LINK_LIBRARIES ${PROJECT_NAME}::Mesh)
# cmake-format: on
