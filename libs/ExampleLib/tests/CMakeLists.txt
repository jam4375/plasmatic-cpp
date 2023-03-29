# cmake-format: off
configure_test_executable(NAME ExampleLibTest
                          SOURCE_FILES main.cpp
                          SOURCE_DIR "."
                          BUILD_LINK_LIBRARIES ${PROJECT_NAME}::ExampleLib)
# cmake-format: on
