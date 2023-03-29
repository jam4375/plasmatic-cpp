include(FetchContent)
FetchContent_Declare(
    cxxopts-external
    GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
    GIT_TAG v2.2.1)

# Don't build unused targets:
set(CXXOPTS_BUILD_EXAMPLES
    OFF
    CACHE BOOL "Set to ON to build examples" FORCE)
set(CXXOPTS_BUILD_TESTS
    OFF
    CACHE BOOL "Set to ON to build tests" FORCE)

FetchContent_MakeAvailable(cxxopts-external)

# Disable warnings originating from external code:
get_target_property(cxxopts_include_dirs cxxopts INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(cxxopts PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${cxxopts_include_dirs}")
