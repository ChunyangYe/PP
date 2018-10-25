cmake_minimum_required(VERSION 3.0)


project(PP_proj)
SET(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "-std=c++11")
if(NOT DEFINED ENV{PARDISO_LIC_PATH})
    message(FATAL_ERROR "not defined environment variable:PARDISO_LIC_PATH")  
endif()
set(SRC_PATH ${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/Include)

add_library(mylib 
			${SRC_PATH}/Parafun.cpp 
			${SRC_PATH}/MeshDefinition.cpp 
			${SRC_PATH}/PardisoSolver.cpp)

add_executable(PP ${SRC_PATH}/main_cmd.cpp)

target_link_libraries(PP 
						mylib
						${PROJECT_SOURCE_DIR}/Library/libOpenMeshTools.a
						${PROJECT_SOURCE_DIR}/Library/libOpenMeshCore.a
						$ENV{PARDISO_LIC_PATH}/libpardiso600-WIN-X86-64.dll)
