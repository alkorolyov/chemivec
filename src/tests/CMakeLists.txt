cmake_minimum_required(VERSION 3.15)
project(chemivec-tests)

# Find Python, Numpy, OpenMP
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_link_libraries(${PROJECT_NAME} chemivec-core)
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
