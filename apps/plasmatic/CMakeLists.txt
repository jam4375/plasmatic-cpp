# Set the version string inside the main source file:
configure_file(main.cpp main.cpp)

# cmake-format: off
configure_executable(NAME plasmatic
                     SOURCE_FILES ${CMAKE_CURRENT_BINARY_DIR}/main.cpp
                     SOURCE_DIR "."
                     BUILD_LINK_LIBRARIES ${PROJECT_NAME}::ProblemTypes cxxopts nlohmann_json::nlohmann_json)
# cmake-format: on

add_test(NAME plasmatic_test COMMAND $<TARGET_FILE:plasmatic> -h)
add_test(NAME plasmatic_test_thermal COMMAND $<TARGET_FILE:plasmatic> -i ${CMAKE_CURRENT_SOURCE_DIR}/config/thermal.json WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
add_test(NAME plasmatic_test_mechanical COMMAND $<TARGET_FILE:plasmatic> -i ${CMAKE_CURRENT_SOURCE_DIR}/config/mechanical.json WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
