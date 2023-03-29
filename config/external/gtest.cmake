if(${PROJECT_NAME}_ENABLE_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG release-1.10.0)

    FetchContent_MakeAvailable(googletest)

    # Disable warnings originating from external code:
    target_compile_options(gtest PRIVATE "-w")
    get_target_property(gtest_include_dirs gtest INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(gtest PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${gtest_include_dirs}")

    target_compile_options(gmock PRIVATE "-w")
    get_target_property(gmock_include_dirs gmock INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(gmock PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${gmock_include_dirs}")

    target_compile_options(gmock_main PRIVATE "-w")
    get_target_property(gmock_main_include_dirs gmock_main INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(gmock_main PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${gmock_main_include_dirs}")

    target_compile_options(gtest_main PRIVATE "-w")
    get_target_property(gtest_main_include_dirs gtest_main INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(gtest_main PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${gtest_main_include_dirs}")

    # Make CMake macros like 'gtest_discover_tests' available:
    include(GoogleTest)
endif()
