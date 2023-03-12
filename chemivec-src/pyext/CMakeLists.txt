cmake_minimum_required(VERSION 3.18)
project(chemivec LANGUAGES C CXX)

# skbuild
#set(CMAKE_MODULE_PATH ${Python3_SITELIB}/skbuild/resources/cmake)
#find_package(PythonExtensions REQUIRED)

# Find Python, Numpy, OpenMP
find_package(Python3 REQUIRED Interpreter Development NumPy)

# Set compiler flags
string(APPEND CMAKE_C_FLAGS " -static")
string(APPEND CMAKE_CXX_FLAGS " -static")
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# Main target, copy pyd ot root dir
python3_add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
set_target_properties (${PROJECT_NAME} PROPERTIES
        LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
        )

message(STATUS ${PROJECT_NAME} : CMAKE_C_FLAGS=${CMAKE_C_FLAGS})
message(STATUS ${PROJECT_NAME} : CMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS})

# PRIVATE because python3_add_library already uses ${PROJECT_NAME} for target_link_libraries
target_link_libraries(${PROJECT_NAME} PRIVATE chemivec-core)