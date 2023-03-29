find_package(Doxygen)

option(${PROJECT_NAME}_BUILD_DOCS "Generate the HTML documentation (requires doxygen)" ${DOXYGEN_FOUND})

if(${PROJECT_NAME}_BUILD_DOCS)
    if(NOT DOXYGEN_FOUND)
        message(FATAL_ERROR "Building documentation requires Doxygen!")
    endif()

    configure_file(Doxyfile.in Doxyfile)

    add_custom_target(
        ${PROJECT_NAME}_docs
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/docs/
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)
endif()
