cmake_minimum_required(VERSION 3.16)
project(chemivec-core LANGUAGES C CXX)


find_package(Python3 COMPONENTS Interpreter Development.Module NumPy REQUIRED)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/src/*.c)

# Main target
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_SOURCES})

target_include_directories(${PROJECT_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src
        PUBLIC ${Python3_INCLUDE_DIRS}
        PUBLIC ${Python3_NumPy_INCLUDE_DIRS}

#        PUBLIC ${Python_INCLUDE_DIRS}
#	    PUBLIC ${Python3_INCLUDE_DIR}
#        PUBLIC ${Python_NumPy_INCLUDE_DIRS}

#        PRIVATE ${OpenMP_C_INCLUDE_DIRS}
#        PRIVATE ${OpenMP_CXX_INCLUDE_DIRS}
        )

target_link_libraries(
        ${PROJECT_NAME}
        PUBLIC indigo-static
        PUBLIC ${Python3_LIBRARIES}
#        PUBLIC ${NumPy_LIBRARIES}
#        PRIVATE OpenMP::OpenMP_C
#        PRIVATE OpenMP::OpenMP_CXX

)