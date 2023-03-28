cmake_minimum_required(VERSION 3.16)

project(chemivec-src LANGUAGES C CXX)

include(FetchContent)
FetchContent_Declare(
        Indigo
        GIT_REPOSITORY https://github.com/epam/Indigo.git
)
FetchContent_MakeAvailable(Indigo)

# Compile flags
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

find_package(OpenMP REQUIRED)
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")

if (UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    string(APPEND CMAKE_C_FLAGS " -static-libstdc++ -static-libgcc")
    string(APPEND CMAKE_CXX_FLAGS " -static-libstdc++ -static-libgcc")
endif()

if (MINGW AND WIN32)
    string(APPEND CMAKE_CXX_FLAGS " -static -Wl,-s")
endif()

message(STATUS CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
message(STATUS CMAKE_C_FLAGS=${CMAKE_C_FLAGS})
message(STATUS CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})
message(STATUS CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
message(STATUS CMAKE_C_FLAGS_RELEASE=${CMAKE_C_FLAGS_RELEASE})
message(STATUS CMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS CIBUILDWHEEL=$ENV{CIBUILDWHEEL})
message(STATUS CMAKE_GENERATOR=${CMAKE_GENERATOR})

add_subdirectory(core)
add_subdirectory(pyext)

if(TEST_CHEMIVEC)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()

