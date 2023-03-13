cmake_minimum_required(VERSION 3.18)
project(chemivec LANGUAGES C CXX)

# Find Python, Numpy
find_package(Python3 REQUIRED Interpreter Development NumPy)

# skbuild
set(CMAKE_MODULE_PATH ${Python3_SITELIB}/skbuild/resources/cmake)
find_package(PythonExtensions REQUIRED)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# Main target, copy pyd to root dir
if(SKBUILD)
    add_library(${PROJECT_NAME} MODULE ${${PROJECT_NAME}_SOURCES})
    python_extension_module(${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME} chemivec-core)
    install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PROJECT_NAME}> ${CMAKE_SOURCE_DIR}
            )
else()
    python3_add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
    set_target_properties (${PROJECT_NAME} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
            )
    target_link_libraries(${PROJECT_NAME} PRIVATE chemivec-core)
endif()






