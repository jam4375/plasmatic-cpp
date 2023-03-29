function(configure_library)
    set(options HEADER_ONLY)
    set(oneValueArgs NAME SOURCE_DIR INTERFACE_DIR)
    set(multiValueArgs SOURCE_FILES BUILD_LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
    cmake_parse_arguments(LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(LibraryTargetName ${PROJECT_NAME}_${LIB_NAME})

    message(
        STATUS
            "Configuring library ${LIB_NAME} with target \"${LibraryTargetName}\" and alias \"${PROJECT_NAME}::${LIB_NAME}\""
    )

    if(${LIB_HEADER_ONLY})
        add_library(${LibraryTargetName} INTERFACE)
        add_library(${PROJECT_NAME}::${LIB_NAME} ALIAS ${LibraryTargetName})

        target_include_directories(${LibraryTargetName} INTERFACE ${LIB_INTERFACE_DIR})

        target_link_libraries(${LibraryTargetName} INTERFACE ${LIB_INTERFACE_LINK_LIBRARIES})
    else()
        add_library(${LibraryTargetName} ${LIB_SOURCE_FILES})
        add_library(${PROJECT_NAME}::${LIB_NAME} ALIAS ${LibraryTargetName})

        target_include_directories(${LibraryTargetName} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/${LIB_SOURCE_DIR})

        target_include_directories(
            ${LibraryTargetName} INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/${LIB_INTERFACE_DIR}>
                                           $<INSTALL_INTERFACE:${LIB_INTERFACE_DIR}> # <prefix>/${LIB_INTERFACE_DIR}
        )

        target_link_libraries(${LibraryTargetName} PRIVATE ${LIB_BUILD_LINK_LIBRARIES})

        target_link_libraries(${LibraryTargetName} PUBLIC ${LIB_INTERFACE_LINK_LIBRARIES})
    endif()
endfunction()
