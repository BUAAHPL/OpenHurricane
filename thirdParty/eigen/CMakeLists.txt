# Author: Rao Sihang
cmake_minimum_required (VERSION 3.18)
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 
 
set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${CMAKE_CURRENT_SOURCE_DIR}/include")