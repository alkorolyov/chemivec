cmake_minimum_required(VERSION 3.18)

project(chemivec-srt LANGUAGES C CXX)

if(POLICY CMP0054)
    # https://cmake.org/cmake/help/v3.9/policy/CMP0054.html
    cmake_policy(SET CMP0054 NEW)
endif()

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

# Set compiler flags
find_package(OpenMP REQUIRED)

if (MINGW)
    # -fwhole breaks the linking
    string(APPEND CMAKE_C_FLAGS " -static ${OpenMP_C_FLAGS} $ENV{CFLAGS}")
    string(APPEND CMAKE_CXX_FLAGS " -static ${OpenMP_CXX_FLAGS} $ENV{CXXFLAGS}")

    set(CMAKE_C_FLAGS_RELEASE " -O3 -ffast-math -DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE " -O3 -ffast-math -DNDEBUG")
elseif(MSVC)
    string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS} $ENV{CFLAGS}")
    string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS} $ENV{CFLAGS}")
    string(APPEND CMAKE_C_FLAGS " -MP -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")
    string(APPEND CMAKE_CXX_FLAGS " -MP -EHs -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")

    set(CMAKE_C_FLAGS_RELEASE " /O2 /fp:fast /DNDEBUG")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /fp:fast /DNDEBUG")
endif()

if (CMAKE_BUILD_TYPE STREQUAL Release)
    string(APPEND CMAKE_C_FLAGS ${CMAKE_C_FLAGS_RELEASE})
    string(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
endif()

message(STATUS CMAKE_C_FLAGS=${CMAKE_C_FLAGS})
message(STATUS CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})


add_subdirectory(core)
add_subdirectory(pyext)

if(ENABLE_TESTS)
    add_subdirectory(tests)
endif()