cmake_minimum_required(VERSION 3.16)
project(chemivec-tests)

enable_testing()

include(CTest)
include(FetchContent)
FetchContent_Declare(
        unity
        GIT_REPOSITORY https://github.com/ThrowTheSwitch/Unity.git
)
FetchContent_MakeAvailable(unity)

file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_link_libraries(${PROJECT_NAME} chemivec-core unity)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
