# Set the version string inside the main source file:
configure_file(main.cpp main.cpp)

# cmake-format: off
configure_executable(NAME plasmatic
                     SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/main.cpp
                     SOURCE_DIR "."
                     BUILD_LINK_LIBRARIES ${PROJECT_NAME}::Mesh cxxopts)
# cmake-format: on

add_test(NAME ${PROJECT_NAME}_plasmatic_test COMMAND $<TARGET_FILE:${PROJECT_NAME}_plasmatic> -h)
