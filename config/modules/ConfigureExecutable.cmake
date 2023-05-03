function(configure_executable)
    set(oneValueArgs NAME SOURCE_DIR)
    set(multiValueArgs SOURCE_FILES BUILD_LINK_LIBRARIES)
    cmake_parse_arguments(EXE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(ExecutableTargetName ${EXE_NAME})

    message(
        STATUS
            "Configuring executable ${EXE_NAME} with target \"${ExecutableTargetName}\" and alias \"${PROJECT_NAME}::${EXE_NAME}\""
    )

    add_executable(${ExecutableTargetName} ${EXE_SOURCE_FILES})
    add_executable(${PROJECT_NAME}::${EXE_NAME} ALIAS ${ExecutableTargetName})

    target_include_directories(${ExecutableTargetName} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${EXE_SOURCE_DIR})

    target_link_libraries(${ExecutableTargetName} PRIVATE ${EXE_BUILD_LINK_LIBRARIES})
endfunction()
