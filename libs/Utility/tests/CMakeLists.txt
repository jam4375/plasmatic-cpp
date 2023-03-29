# cmake-format: off
configure_test_executable(NAME UtilityTest
                          SOURCE_FILES main.cpp
                          SOURCE_DIR "."
                          BUILD_LINK_LIBRARIES ${PROJECT_NAME}::Utility)
# cmake-format: on
