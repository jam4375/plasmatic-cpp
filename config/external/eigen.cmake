include(FetchContent)
FetchContent_Declare(
  Eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE
  GIT_PROGRESS TRUE)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "Set to ON to build docs" FORCE)
set(EIGEN_BUILD_PKGCONFIG OFF)
set(BUILD_TESTING OFF CACHE BOOL "Set to ON to build tests" FORCE)

FetchContent_MakeAvailable(Eigen)

get_target_property(EIGEN_INCLUDES eigen INTERFACE_INCLUDE_DIRECTORIES)
set_target_properties(eigen PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${EIGEN_INCLUDES}")