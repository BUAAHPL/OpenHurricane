# Author: Rao Sihang
# CMakeList.txt : CMake project for OpenHurricane, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project (
         chemSolvers
         LANGUAGES CXX
        )
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 

#get_code_files(${CMAKE_CURRENT_SOURCE_DIR})

#get_property(local_prop GLOBAL PROPERTY source_list) 
get_property(tmp_src_list GLOBAL PROPERTY global_src_list) 
get_property(tmp_header_list GLOBAL PROPERTY global_header_list) 
get_property(tmp_include_dirs GLOBAL PROPERTY global_include_dirs1)

set(src_files "${CMAKE_CURRENT_SOURCE_DIR}/chemSolvers.cpp")
list ( APPEND tmp_src_list  ${src_files} )     

message(STATUS "Setting project: ${PROJECT_NAME}")
# Add source to this project's executable.
add_executable (${PROJECT_NAME} ${tmp_src_list} ${tmp_header_list})

get_property(PRJ_tmp_LIBRARIES GLOBAL PROPERTY PRJ_LIBRARIES) 
get_property(PRJ_tmp_COMPILE_FEATURES GLOBAL PROPERTY PRJ_COMPILE_FEATURES) 
get_property(PRJ_tmp_COMPILE_DEF GLOBAL PROPERTY PRJ_COMPILE_DEF) 


target_link_libraries( ${PROJECT_NAME}
	PRIVATE
		${PRJ_tmp_LIBRARIES}
)

target_include_directories( ${PROJECT_NAME} 
	PRIVATE 
		${tmp_include_dirs}
)

target_compile_features ( ${PROJECT_NAME} 
	PRIVATE 
		${PRJ_tmp_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
	PRIVATE
	   ${PRJ_tmp_COMPILE_DEF} 
)
#message(STATUS "${PROJECT_NAME} using libraries: " ${PRJ_tmp_LIBRARIES})
#message(STATUS "Compile definitions of ${PROJECT_NAME}: ${PRJ_tmp_COMPILE_DEF} ")
# install 
install(TARGETS ${PROJECT_NAME} DESTINATION bin)

