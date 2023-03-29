# cmake-format: off
configure_library(NAME ExampleLib
                  SOURCE_FILES ExampleLib.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES ""
                  INTERFACE_LINK_LIBRARIES "")
# cmake-format: on

add_subdirectory(tests)
