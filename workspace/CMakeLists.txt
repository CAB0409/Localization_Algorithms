cmake_minimum_required (VERSION 3.0)

project(tdoa)
cmake_policy(SET CMP0054 OLD)

#find packages
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
find_package(Eigen3 3.3 REQUIRED)
#find_package(ZMQ REQUIRED)
find_package(Armadillo REQUIRED)
find_package(Matlab REQUIRED)

#Compiler Definitions
add_definitions(-std=c++11)

#subdirectories
#message("source directory: " + ${PROJECT_SOURCE_DIR})

###
# Group Files Together
# These will be linked to DARC exe
###
#Group all lib files
file(GLOB_RECURSE PROJECT_SRCS ${PROJECT_SOURCE_DIR}/lib/*.cpp)

#include_directories(${ZMQ_INCLUDE_DIR})
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/lib)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${Matlab_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${PROJECT_SRCS})

###
# Library Dependencies
###

#Proto dependencies
if(PROTO_FILES)
    add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_protolib)
    target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_protolib)
endif(PROTO_FILES)
#Waveform dependencies
if(WAVEFORM_SRCS)
    add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_waveformlib)
    target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_waveformlib)
endif(WAVEFORM_SRCS)

target_link_libraries(${PROJECT_NAME}
                      ${ARMADILLO_LIBRARY}
                      ${Matlab_MEX_LIBRARY}
                      Eigen3::Eigen)
