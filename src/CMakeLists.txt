cmake_minimum_required(VERSION 3.16)

project(chemivec-src LANGUAGES C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()



# Compile flags
find_package(OpenMP REQUIRED)
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS} $ENV{CFLAGS}")
string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS} $ENV{CXXFLAGS}")

if (UNIX OR MINGW)
    string(APPEND CMAKE_C_FLAGS " -fvisibility=hidden")
    string(APPEND CMAKE_CXX_FLAGS " -fvisibility=hidden -fvisibility-inlines-hidden")

    if (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
        string(APPEND CMAKE_CXX_FLAGS " -static-libstdc++")
        if (MINGW)
            string(APPEND CMAKE_CXX_FLAGS " -static-libgcc")
        endif()
        if (WIN32)
            string(APPEND CMAKE_CXX_FLAGS " -static")
        endif()
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        string(APPEND CMAKE_CXX_FLAGS " -Wall -Wextra -stdlib=libc++")
    endif()

    if(MSYS OR CYGWIN OR MINGW)
        set(CMAKE_SHARED_LIBRARY_PREFIX "")
        if (MSYS OR CYGWIN)
            string(APPEND CMAKE_C_FLAGS " -D_GNU_SOURCE")
            string(APPEND CMAKE_CXX_FLAGS " -D_GNU_SOURCE")
        endif()
    endif()
elseif (MSVC)
    string(APPEND CMAKE_C_FLAGS " -MP -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")
    string(APPEND CMAKE_CXX_FLAGS " -MP -EHs -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")
endif ()


#if (MINGW)
#    # -fwhole breaks the linking
#    string(APPEND CMAKE_C_FLAGS " -static $ENV{CFLAGS}")
#    string(APPEND CMAKE_CXX_FLAGS " -static $ENV{CXXFLAGS}")
#
#    set(CMAKE_C_FLAGS_RELEASE " -O3 -ffast-math -DNDEBUG")
#    set(CMAKE_CXX_FLAGS_RELEASE " -O3 -ffast-math -DNDEBUG")
#elseif(MSVC)
#    string(APPEND CMAKE_C_FLAGS " -MP -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")
#    string(APPEND CMAKE_CXX_FLAGS " -MP -EHs -D_CRT_SECURE_NO_WARNINGS -D_WIN32_WINNT=0x0601 -DWINVER=0x0601")
#
#    set(CMAKE_C_FLAGS_RELEASE " /O2 /fp:fast /DNDEBUG")
#    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /fp:fast /DNDEBUG")
#endif()

#if (CMAKE_BUILD_TYPE STREQUAL Release)
#    string(APPEND CMAKE_C_FLAGS ${CMAKE_C_FLAGS_RELEASE})
#    string(APPEND CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE})
#endif()

message(STATUS CMAKE_C_FLAGS=${CMAKE_C_FLAGS})
message(STATUS CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})
message(STATUS CIBUILDWHEEL=${CIBUILDWHEEL})
message(STATUS CMAKE_GENERATOR=${CMAKE_GENERATOR})

if(CMAKE_GENERATOR STREQUAL "Unix Makefiles")
    set(CMAKE_GENERATOR Ninja)
endif()

add_subdirectory(core)
add_subdirectory(pyext)

if(TEST_CHEMIVEC)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif()