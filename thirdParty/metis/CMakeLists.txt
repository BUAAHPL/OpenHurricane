# Author: Rao Sihang
cmake_minimum_required (VERSION 3.18)
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" )

if(NOT WIN32)# Linux
    if( CMAKE_BUILD_TYPE STREQUAL "Debug" )        
        get_filename_component(lib_abs "lib/Linux/libmetisd.a" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
    else( CMAKE_BUILD_TYPE STREQUAL "Debug" )#release
        get_filename_component(lib_abs "lib/Linux/libmetis.a" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
    endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
     
else(NOT WIN32)#Windows
    if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
        get_filename_component(lib_abs "lib/Windows/metisd.lib" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
    else( CMAKE_BUILD_TYPE STREQUAL "Debug" )#release
        get_filename_component(lib_abs "lib/Windows/metis.lib" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
    endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )

endif(NOT WIN32)
set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${CMAKE_CURRENT_SOURCE_DIR}/include")
