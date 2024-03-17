# Author: Rao Sihang
cmake_minimum_required (VERSION 3.18)
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 

if(NOT WIN32)# Linux
        # If the program uses CGNS, add the libraries of CGNS to the program.
    if(CGNS_ENABLE)    
        get_filename_component(lib_abs "lib/Linux/lib/libcgns.a" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
        set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${CMAKE_CURRENT_SOURCE_DIR}/include")
    endif(CGNS_ENABLE)

else(NOT WIN32)#Windows
    # If the program uses CGNS, add the libraries of CGNS to the program.
    if(CGNS_ENABLE)    
        get_filename_component(lib_abs "lib/Windows/lib/cgns.lib" ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${lib_abs}")
        set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${CMAKE_CURRENT_SOURCE_DIR}/include")
    endif(CGNS_ENABLE)
endif(NOT WIN32)