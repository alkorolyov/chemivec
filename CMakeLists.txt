cmake_minimum_required(VERSION 3.15)

project(Chemivec LANGUAGES NONE)


# Override Indigo options to install indigo-core only
set(BUILD_INDIGO_WRAPPERS OFF CACHE STRING "Build all existing Indigo API wrappers" FORCE)
set(BUILD_INDIGO_UTILS OFF CACHE STRING "Build indigo-cano, indigo-deco and indigo-depict utilities" FORCE)
set(BUILD_BINGO OFF CACHE STRING "Build all existing Bingo cartridges" FORCE)
set(BUILD_BINGO_ELASTIC OFF CACHE STRING "Build Bingo Elastic cartridges for Java and Python" FORCE)
set(WITH_STATIC ON CACHE STRING "Build Indigo static library as well as shared" FORCE)
set(ENABLE_TESTS OFF CACHE STRING "Enable Indigo CTest tests" FORCE)

# Tests

set(TEST_CHEMIVEC ON CACHE STRING "Enable Chemivec CTest tests")
add_subdirectory(src)
