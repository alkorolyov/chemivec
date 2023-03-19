cmake_minimum_required(VERSION 3.16)
project(chemivec-core LANGUAGES C CXX)

# Download Indigo package from git
include(FetchContent)
FetchContent_Declare(
        Indigo
        GIT_REPOSITORY https://github.com/epam/Indigo.git
)
FetchContent_MakeAvailable(Indigo)

# Python, Numpy, OpenMP



if(UNIX)
    execute_process(COMMAND yum install -y python3-devel python3-numpy)

find_package(Python3 REQUIRED Interpreter Development NumPy)
find_package(OpenMP REQUIRED)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

# Main target
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})

target_include_directories(${PROJECT_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
        PUBLIC ${Python3_INCLUDE_DIRS}
        PUBLIC ${Python3_NumPy_INCLUDE_DIRS}
        PUBLIC ${OpenMP_C_INCLUDE_DIRS}
        PUBLIC ${OpenMP_CXX_INCLUDE_DIRS}
        )


target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC indigo-static
        PUBLIC ${Python3_LIBRARIES}
        PUBLIC ${NumPy_LIBRARIES}
)