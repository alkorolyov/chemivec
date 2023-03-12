cmake_minimum_required(VERSION 3.18)
project(chemivec-src LANGUAGES NONE)

# Override Indigo options to install indigo-core only
set(BUILD_INDIGO_WRAPPERS OFF CACHE STRING "Build all existing Indigo API wrappers" FORCE)
set(BUILD_INDIGO_UTILS OFF CACHE STRING "Build indigo-cano, indigo-deco and indigo-depict utilities" FORCE)
set(BUILD_BINGO OFF CACHE STRING "Build all existing Bingo cartridges" FORCE)
set(BUILD_BINGO_ELASTIC OFF CACHE STRING "Build Bingo Elastic cartridges for Java and Python" FORCE)
set(WITH_STATIC ON CACHE STRING "Build Indigo static library as well as shared" FORCE)

# Tests
set(ENABLE_TESTS ON CACHE STRING "Enable CTest tests" FORCE)

add_subdirectory(chemivec-src)
