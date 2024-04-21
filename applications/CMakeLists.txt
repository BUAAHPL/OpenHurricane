# ---------------------------------------------------------------
# Programmer(s): Rao Sihang
# ---------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)

message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 
# The main application of Hurricane must be build.
add_subdirectory (HurricaneMain)

option(APPLICATIONS_chemSolvers_ENABLE "If should build application: chemSolvers " OFF)
if(APPLICATIONS_chemSolvers_ENABLE)
	add_subdirectory (chemSolvers)
endif(APPLICATIONS_chemSolvers_ENABLE)