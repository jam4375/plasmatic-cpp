FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz)
FetchContent_MakeAvailable(json)

# Disable warnings originating from external code:
get_target_property(JSON_INCLUDES nlohmann_json INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(nlohmann_json PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${JSON_INCLUDES}")