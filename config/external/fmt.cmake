FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG master
)
FetchContent_MakeAvailable(fmt)

# Disable warnings originating from external code:
get_target_property(FMT_INCLUDES fmt INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(fmt PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${FMT_INCLUDES}")
set_target_properties(fmt PROPERTIES EXCLUDE_FROM_ALL ON)

get_target_property(FMT_HEADER_ONLY_INCLUDES fmt-header-only INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(fmt-header-only PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${FMT_HEADER_ONLY_INCLUDES}")