cmake_minimum_required(VERSION 3.18)
project(chemivec-tests)

# Find Python, Numpy, OpenMP
string(APPEND CMAKE_C_FLAGS " ${OpenMP_C_FLAGS}")
string(APPEND CMAKE_CXX_FLAGS " ${OpenMP_CXX_FLAGS}")

file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

add_executable(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})

#target_include_directories(${PROJECT_NAME} PUBLIC
#        ${Python3_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR}
#        )

target_link_libraries(${PROJECT_NAME} chemivec-core)
