cmake_minimum_required(VERSION 2.8)
 
PROJECT(aerosol_in_stratosphere)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(aerosol_in_stratosphere aerosol_in_stratosphere.cpp)
 
target_link_libraries(aerosol_in_stratosphere ${VTK_LIBRARIES})
