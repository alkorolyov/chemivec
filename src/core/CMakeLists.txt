cmake_minimum_required(VERSION 3.16)
project(chemivec-core LANGUAGES C CXX)

# Download Indigo package from git
include(FetchContent)
FetchContent_Declare(
        Indigo
        GIT_REPOSITORY https://github.com/epam/Indigo.git
)
FetchContent_MakeAvailable(Indigo)

include(FindPython3)
if(WIN32)
    find_package(Python3 REQUIRED Interpreter Development NumPy)
#    set(Python_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
#    set(Python_LIBRARIES ${Python3_LIBRARIES})
#    set(Python_NumPy_INCLUDE_DIRS ${Python3_NumPy_INCLUDE_DIRS})
elseif(UNIX AND SKBUILD)
    find_package(Python3 REQUIRED)
endif()
find_package(OpenMP REQUIRED)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

# Main target
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})

target_include_directories(${PROJECT_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
        PUBLIC ${Python_INCLUDE_DIRS}
        PUBLIC ${Python3_INCLUDE_DIRS}
	PUBLIC ${Python3_INCLUDE_DIR}
        PUBLIC ${Python_NumPy_INCLUDE_DIRS}
        PUBLIC ${Python3_NumPy_INCLUDE_DIRS}
#        PUBLIC ${OpenMP_C_INCLUDE_DIRS}
#        PUBLIC ${OpenMP_CXX_INCLUDE_DIRS}
        )

target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC indigo-static
        PUBLIC ${Python_LIBRARIES}
        PUBLIC ${Python3_LIBRARIES}
#        PUBLIC ${NumPy_LIBRARIES}
)