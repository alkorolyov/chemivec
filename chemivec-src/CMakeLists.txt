cmake_minimum_required(VERSION 3.18)

add_subdirectory(core)
add_subdirectory(pyext)

if(ENABLE_TESTS)
    add_subdirectory(tests)
endif()