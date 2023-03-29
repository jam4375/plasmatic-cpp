cmake_minimum_required(VERSION 3.20.0)

project(plasmatic CXX)

# Setup project configuration:
include("${${PROJECT_NAME}_SOURCE_DIR}/config/config.cmake")

# Add documentation directory:
add_subdirectory(docs)

# Add source code directories:
add_subdirectory(apps)
add_subdirectory(libs)
