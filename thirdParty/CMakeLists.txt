# Author: Rao Sihang
cmake_minimum_required (VERSION 3.18)
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 

add_subdirectory(hdf5)
add_subdirectory(metis)
add_subdirectory(eigen)
add_subdirectory(cgns)