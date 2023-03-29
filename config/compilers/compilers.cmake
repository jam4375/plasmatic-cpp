set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(EMSCRIPTEN)
    # Output to html if using EMSCRIPTEN
    set(CMAKE_EXECUTABLE_SUFFIX ".html")

    # Disable clang-tidy when using EMSCRIPTEN
    set(CMAKE_CXX_CLANG_TIDY
        ""
        CACHE INTERNAL "Disable clang-tidy when using EMSCRIPTEN")
endif()

if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(STATUS "Loading Clang cmake settings")
    include("${CMAKE_CURRENT_LIST_DIR}/Clang.cmake")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Loading GCC cmake settings")
    include("${CMAKE_CURRENT_LIST_DIR}/GCC.cmake")
else()
    message(FATAL_ERROR "Unsupported compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()
