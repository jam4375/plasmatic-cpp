# cmake-format: off
configure_library(NAME Utility
                  SOURCE_FILES Utility.cpp ExecutablePath.cpp
                  SOURCE_DIR "."
                  INTERFACE_DIR "interface"
                  BUILD_LINK_LIBRARIES ""
                  INTERFACE_LINK_LIBRARIES "fmt::fmt-header-only")
# cmake-format: on

add_subdirectory(tests)
