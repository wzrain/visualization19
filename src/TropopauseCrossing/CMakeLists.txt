cmake_minimum_required(VERSION 2.8)
 
PROJECT(tropopause_crossing)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(tropopause_crossing tropopause_crossing.cpp)
 
target_link_libraries(tropopause_crossing ${VTK_LIBRARIES})
