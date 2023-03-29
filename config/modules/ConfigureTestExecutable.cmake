function(configure_test_executable)
    if(NOT ${PROJECT_NAME}_ENABLE_TESTS)
        return()
    endif()

    set(oneValueArgs NAME SOURCE_DIR)
    set(multiValueArgs SOURCE_FILES BUILD_LINK_LIBRARIES)
    cmake_parse_arguments(TEST_EXE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(TestExecutableTargetName ${PROJECT_NAME}_${TEST_EXE_NAME})

    message(
        STATUS
            "Configuring test executable ${TEST_EXE_NAME} with target \"${TestExecutableTargetName}\" and alias \"${PROJECT_NAME}::${TEST_EXE_NAME}\""
    )

    add_executable(${TestExecutableTargetName} ${TEST_EXE_SOURCE_FILES})
    add_executable(${PROJECT_NAME}::${TEST_EXE_NAME} ALIAS ${TestExecutableTargetName})

    target_include_directories(${TestExecutableTargetName} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_EXE_SOURCE_DIR})

    target_link_libraries(${TestExecutableTargetName} PRIVATE ${TEST_EXE_BUILD_LINK_LIBRARIES} gtest)

    # disable warning that's fundamental to gtest:
    target_compile_options(${TestExecutableTargetName} PRIVATE "-Wno-global-constructors")

    # disable clang-tidy for the test target (otherwise it will complain about gtest's macros and it doesn't look like
    # they can be cleanly ignored)
    set_target_properties(${TestExecutableTargetName} PROPERTIES CXX_CLANG_TIDY "")

    gtest_discover_tests(${TestExecutableTargetName})
endfunction()
