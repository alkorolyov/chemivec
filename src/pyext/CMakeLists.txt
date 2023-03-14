cmake_minimum_required(VERSION 3.18)
project(_chemivec LANGUAGES C CXX)

# Find Python
find_package(Python3 REQUIRED Development Interpreter NumPy)

# skbuild
set(CMAKE_MODULE_PATH ${Python3_SITELIB}/skbuild/resources/cmake)
find_package(PythonExtensions REQUIRED)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

# Main target
if(SKBUILD)
    add_library(${PROJECT_NAME} MODULE ${${PROJECT_NAME}_SOURCES})
    python_extension_module(${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME} chemivec-core)
    install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION src/chemivec)

else()
    python3_add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
    set_target_properties (${PROJECT_NAME} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}
            )
    target_link_libraries(${PROJECT_NAME} PRIVATE chemivec-core)
endif()





