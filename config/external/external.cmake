# Temporarily disable clang-tidy configuration:
set(CMAKE_CXX_CLANG_TIDY_OLD ${CMAKE_CXX_CLANG_TIDY})
set(CMAKE_CXX_CLANG_TIDY
    ""
    CACHE INTERNAL "Disable clang-tidy for external targets")

include("${CMAKE_CURRENT_LIST_DIR}/gtest.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/cxxopts.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/fmt.cmake")
include("${CMAKE_CURRENT_LIST_DIR}/eigen.cmake")

# Restore clang-tidy configuration:
set(CMAKE_CXX_CLANG_TIDY
    ${CMAKE_CXX_CLANG_TIDY_OLD}
    CACHE BOOL "Re-enable clang-tidy for ${PROJECT_NAME}" FORCE)
