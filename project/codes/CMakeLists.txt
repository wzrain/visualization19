cmake_minimum_required(VERSION 2.8)
 
PROJECT(volcano)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(volcano src/volcano.cpp)
add_executable(volcanoArgv src/volcanoArgv.cpp)
 
target_link_libraries(volcano ${VTK_LIBRARIES})
target_link_libraries(volcanoArgv ${VTK_LIBRARIES})