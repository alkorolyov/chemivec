cmake_minimum_required(VERSION 3.16)
project(_chemivec LANGUAGES C CXX)

# Find all sources
file(GLOB ${PROJECT_NAME}_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/*.c)


## Python extension naming
#add_library(${PROJECT_NAME} MODULE ${${PROJECT_NAME}_SOURCES})
#target_link_libraries(${PROJECT_NAME} PRIVATE chemivec-core)
#set_target_properties(${PROJECT_NAME} PROPERTIES PREFIX "")
#
#if(WIN32)
#    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".pyd")
#else()
#    set_target_properties(${PROJECT_NAME} PROPERTIES SUFFIX ".so")
#endif()

# Find Python
find_package(Python3 REQUIRED Interpreter Development.Module)
Python3_add_library(${PROJECT_NAME} ${${PROJECT_NAME}_SOURCES})
target_link_libraries(${PROJECT_NAME} PRIVATE chemivec-core)

if(SKBUILD)
    install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/src/chemivec)
    set_target_properties (${PROJECT_NAME} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/src/chemivec
            )

# Clion CMake build
else()
    set_target_properties (${PROJECT_NAME} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/src/chemivec
            )
    # Run pytest after build
    if (TEST_CHEMIVEC)
        add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
                WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                COMMAND pytest ./tests -v)
    endif()
endif()







