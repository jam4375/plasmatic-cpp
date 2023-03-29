# Append to module path:
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}")

# Add custom cmake modules:
include(ConfigureLibrary)
include(ConfigureExecutable)
include(ConfigureTestExecutable)
