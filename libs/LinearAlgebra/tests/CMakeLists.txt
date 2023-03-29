# cmake-format: off
configure_test_executable(NAME LinearAlgebraTest
                          SOURCE_FILES main.cpp
                          SOURCE_DIR "."
                          BUILD_LINK_LIBRARIES ${PROJECT_NAME}::LinearAlgebra)
# cmake-format: on