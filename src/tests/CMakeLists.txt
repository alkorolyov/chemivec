cmake_minimum_required(VERSION 3.15)
project(chemivec-tests)


#set_property(GLOBAL PROPERTY CTEST_TARGETS_ADDED 1)
enable_testing()

include(CTest)
include(FetchContent)
FetchContent_Declare(
        unity
        GIT_REPOSITORY https://github.com/ThrowTheSwitch/Unity.git
)
FetchContent_MakeAvailable(unity)


# Find Python, Numpy, OpenMP
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

message(STATUS "Tests added as: ${PROJECT_NAME}.exe")

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_link_libraries(${PROJECT_NAME} chemivec-core unity)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
