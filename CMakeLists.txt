# ---------------------------------------------------------------
# Programmer(s): Rao Sihang
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Top-level CMake project file, do global configuration
# and include sub-projects here.
# ---------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)
include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenHurricane.cmake" )

set_property(GLOBAL PROPERTY PRJ_COMPILE_DEF )
set_property(GLOBAL PROPERTY PRJ_LIBRARIES )
set_property(GLOBAL PROPERTY PRJ_COMPILE_FEATURES )

set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_FEATURES cxx_std_17)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set (project_homepage "http://http://172.16.104.17:9818/summary/OpenHurricane-master-v1.0.0.git")

#======================================================
# 1. Setting the project info:
#    (1) Project name
#    (2) Project version
#    (3) Project description
#======================================================
project (
         OpenHurricane 
         VERSION 2.0.0
         DESCRIPTION "Highly Universal Rocket & Ramjet sImulation Codes for ANalysis and Evaluation"
         LANGUAGES CXX
        )

if (POLICY CMP0074)
  cmake_policy (SET CMP0074 NEW)
endif ()

if (POLICY CMP0083)
  cmake_policy (SET CMP0083 NEW)
endif ()

message ( STATUS "Welcome to build this program. The program information is as follows:")
message ( STATUS "Project name          : ${CMAKE_PROJECT_NAME}"        )
message ( STATUS "Major version         : ${PROJECT_VERSION_MAJOR}"     )
message ( STATUS "Minor version         : ${PROJECT_VERSION_MINOR}"     )
message ( STATUS "Subminor version      : ${PROJECT_VERSION_PATCH}"     )
message ( STATUS "OpenHurricane description : ${CMAKE_PROJECT_DESCRIPTION}" )
message ( STATUS "The cmake file is written by Rao Sihang." )

#======================================================
# 2. Disallow in-source builds.
#    
#======================================================
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed. Please make a new directory (called a build directory) and run CMake from there. You may need to remove CMakeCache.txt. ")
endif()

#======================================================
# 3. Generating configure file of OpenHurricane
#    
#======================================================
message( STATUS "Generating configure file OpenHurricaneConfig.hpp in ${PROJECT_BINARY_DIR}")
configure_file(cmake/OpenHurricaneConfig.h.in OpenHurricaneConfig.hpp)
get_property(tg_include_dirs GLOBAL PROPERTY global_include_dirs1) 
if("${tg_include_dirs}" STREQUAL "")
    set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${PROJECT_BINARY_DIR}")
else("${tg_include_dirs}" STREQUAL "")
    set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${PROJECT_BINARY_DIR}")
endif("${tg_include_dirs}" STREQUAL "")

#======================================================
# 4. Set install dir
#    
#======================================================
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/OpenHurricane CACHE PATH "My default install prefix" FORCE)
    message(STATUS "Setting default install prefix as " ${CMAKE_INSTALL_PREFIX})
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    #set(CMAKE_INSTALL_PREFIX ${CMAKE_BINARY_DIR}/OpenHurricane CACHE PATH "My default install prefix" FORCE)
    message(STATUS "Setting install prefix to " ${CMAKE_INSTALL_PREFIX})
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
message(STATUS "CAMKE_BINARY_DIR is " ${CMAKE_BINARY_DIR})

#======================================================
# 5. Set precision: float or double
#    
#======================================================
set(DOCSTR "Single or double")
hurricane_option(HURRICANE_PRECISION STRING "${DOCSTR}" "DOUBLE")
if(HURRICANE_PRECISION MATCHES "DOUBLE")
    # Double precision
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HURRICANE_DP)
    message(STATUS "The OpenHurricane program is using double precision")
    set(Precision_dir "double")
elseif(HURRICANE_PRECISION MATCHES "SINGLE")
    # Single precision
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HURRICANE_SP)
    message(STATUS "The OpenHurricane program is using single precision")
    set(Precision_dir "float")
else()
    message(FATAL_ERROR "Invalid precision type")
endif()

#======================================================
# 6. Set compiler flags
#    
#======================================================
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
if(NOT WIN32)# Linux
    set(CMAKE_BINARY_DIR "${CMAKE_BINARY_DIR}/Linux")
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wno-unused-variable -pthread -Wclass-conversion -fpermissive -Wunused-function -Wpedantic")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -ldl")
        set(parallel_dir "${parallel_dir}Debug")
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${parallel_dir}/${Precision_dir})
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release" )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pthread -fopenmp -Wclass-conversion -fpermissive -mavx2 -m64")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed -ldl")
        set(parallel_dir "${parallel_dir}Release")
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${parallel_dir}/${Precision_dir})
    endif(CMAKE_BUILD_TYPE STREQUAL "Debug" )
else(NOT WIN32)#Windows
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W3 /permissive- /MP /GS /Zc:wchar_t /Gm- /Od /sdl /wd4819")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Gd /MDd /diagnostics:column ")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  /DEBUG:FASTLINK /MACHINE:X64 /INCREMENTAL /ERRORREPORT:PROMPT /NOLOGO /TLBID:1 ")        
    elseif(CMAKE_BUILD_TYPE STREQUAL "Release" )
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /Gw /GL /Gy /GA /Gd /O2 /Ob2 /Ot /Oi- /MD /MP /wd4819 /fp:fast /Qpar-report:2 /QIntel-jcc-erratum /W3 /FC /sdl /diagnostics:column /Zc:inline -openmp:experimental -Qvec-report:2 /arch:AVX2")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST /DYNAMICBASE /NXCOMPAT /LTCG:incremental /DEBUG:NONE /MACHINE:X64 /INCREMENTAL:NO /ERRORREPORT:PROMPT /OPT:REF,ICF,LBR /NOLOGO /TLBID:1")
    elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  /GS /GL /W3 /Gy /Zc:wchar_t /Zi /Gm- /O2 /wd4819 /sdl -openmp:experimental")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:inline /fp:fast /MP /diagnostics:column /FC /Fa /errorReport:prompt /WX- /Zc:forScope")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Gd /Oi- /MD /FC /Ot")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /MANIFEST /DYNAMICBASE /NXCOMPAT /LTCG:incremental /DEBUG:NONE /MACHINE:X64 /INCREMENTAL:NO /ERRORREPORT:PROMPT /OPT:REF,ICF,LBR /NOLOGO /TLBID:1")
    endif( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    set(CMAKE_BINARY_DIR "${CMAKE_BINARY_DIR}/Windows")        
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${Precision_dir})
endif(NOT WIN32)
if(NOT MSVC)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
        #add_compile_definitions(_DEBUG)
        set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF _DEBUG)
    endif()
else(NOT MSVC)# 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc /wd4127 /wd4505 /wd4714 /permissive-")
    #add_compile_definitions(_CRT_SECURE_NO_DEPRECATE)
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF _CRT_SECURE_NO_DEPRECATE)
endif(NOT MSVC)

if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
    #add_compile_definitions(_DEBUG)
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HUR_DEBUG)
    option(FULL_DEBUG_ENABLE "If should build with full-debug" OFF)
    if(FULL_DEBUG_ENABLE)
        set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF FULLDEBUG)
    endif(FULL_DEBUG_ENABLE)
endif()

#======================================================
# 7. Set EIGEN
#If you want to guarantee that the Eigen code that you are #including is licensed
#under the MPL2 and possibly more permissive licenses (like BSD), #define this
#preprocessor symbol:
#   EIGEN_MPL2_ONLY
#======================================================
#add_compile_definitions(EIGEN_MPL2_ONLY)
set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF EIGEN_MPL2_ONLY)

#======================================================
# 8. Set MPI
#    
#======================================================
option(MPI_ENABLE "If the program uses MPI" ON)
if(MPI_ENABLE)
	# Have MPI
    message(STATUS "Finding MPI...")
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
	    #add_compile_definitions(MPI_PARALLEL)
        set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF MPI_PARALLEL)
	    set(parallel_dir "parallel")
        message(STATUS "The MPI will be used")
        include_directories(${MPI_INCLUDE_PATH})
        message(STATUS "MPI include path: " ${MPI_INCLUDE_PATH})        
        message(STATUS "MPI libraries path: " ${MPI_LIBRARIES})
        set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${MPI_LIBRARIES}")
    else(MPI_FOUND)        
        message(STATUS "MPI not found. The MPI would not be used.")
    endif(MPI_FOUND)
endif(MPI_ENABLE)

#======================================================
# 9. Set CVODE
#    
#======================================================
#option(CVODE_ENABLE "If the program uses CVODE" ON)
#if(CVODE_ENABLE)
#    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF USES_CVODE)
#    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF SUNDIALS_STATIC_DEFINE)
#    message(STATUS "CVODE will be compiled into this program.")
#endif(CVODE_ENABLE)

#======================================================
# 10. Set CGNS
#    
#======================================================
option(CGNS_ENABLE "If the program uses CGNS" ON)
if(CGNS_ENABLE)
    #add_compile_definitions(USES_CGNS)
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF USES_CGNS)
    message(STATUS "CGNS will be compiled into this program.")
endif(CGNS_ENABLE)

#======================================================
# 11. Set CUDA
#    
#======================================================
option(CUDA_ENABLE "If the program uses CUDA" OFF)
if(CUDA_ENABLE)
	# Have CUDA
    message(STATUS "Finding CUDA...")
    include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/OpenHurricaneSetupCUDA.cmake" )
endif(CUDA_ENABLE)

#======================================================
# 12. Setup integer size
#    
#======================================================
set(DOCSTR "Signed 64-bit (64) or signed 32-bit (32) integer")
hurricane_option(HURRICANE_INDEX_SIZE STRING "${DOCSTR}" "32")
if(HURRICANE_INDEX_SIZE MATCHES "32")
    # Double precision
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HURRICANE_32_INT)
    message(STATUS "The OpenHurricane program is using 32-bit (32) integer")
elseif(HURRICANE_INDEX_SIZE MATCHES "64")
    # Single precision
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HURRICANE_64_INT)
    message(STATUS "The OpenHurricane program is using 64-bit (64) integer")
else()
    message(FATAL_ERROR "Invalid index size")
endif()

#======================================================
# 13. Test process time
#    
#======================================================
option(TEST_PROCESS_TIME_ENABLE "Whether to test the process time of the program " OFF)
if(TEST_PROCESS_TIME_ENABLE)
    message(STATUS "Setting the TEST_PROCESS_TIME_ENABLE flag to ON. And the process time of the OpenHurricane would be tested")
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF TEST_PROCESS_TIME)  
    option(TEST_PROCESS_TIME_ENABLE_DEBUG "Whether to debug the process time of the program " OFF)
    if(TEST_PROCESS_TIME_ENABLE_DEBUG)
        message(STATUS "Setting the TEST_PROCESS_TIME_ENABLE_DEBUG flag to ON. And the debug process time of the OpenHurricane would be tested")
        set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF TEST_PROCESS_TIME_DEBUG)    
    endif(TEST_PROCESS_TIME_ENABLE_DEBUG)  
endif(TEST_PROCESS_TIME_ENABLE)

message ( STATUS "Welcome to build project ${PROJECT_NAME}." )
message ( STATUS "The Operator-System is ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}." )
message ( STATUS "The CMAKE_SOURCE_DIR is ${CMAKE_SOURCE_DIR}." )
message ( STATUS "The CMAKE_BINARY_DIR is ${CMAKE_BINARY_DIR}." ) 
message ( STATUS "The working directory is ${CMAKE_BINARY_DIR}." )
message ( STATUS "The cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 
message ( STATUS "This current binary directory is ${CMAKE_CURRENT_BINARY_DIR}." )
message ( STATUS "This current source directory is ${CMAKE_CURRENT_SOURCE_DIR}." )
message ( STATUS "The C compiler is ${CMAKE_C_COMPILER_ID}" )
message ( STATUS "The CXX compiler is ${CMAKE_CXX_COMPILER_ID}" )
message ( STATUS "The CMAKE_GENERATOR is ${CMAKE_GENERATOR}" )
message ( STATUS "The CMAKE_CXX_COMPILER_VERSION is ${CMAKE_CXX_COMPILER_VERSION}" )
message ( STATUS "Build Type: " ${CMAKE_BUILD_TYPE} ${CMAKE_CXX_FLAGS})
message ( STATUS "Link flag: " ${CMAKE_EXE_LINKER_FLAGS})
message ( STATUS "The executable binary file will be written in: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")


#======================================================
# 14. Logger option: None, Less or Full
#======================================================
set(DOCSTR " None, Less or Full")
hurricane_option(HURRICANE_LOGGER STRING "${DOCSTR}" "LESS")
if(HURRICANE_LOGGER MATCHES "NONE")
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HUR_NONE_LOGGER)
    message(STATUS "The OpenHurricane program is using none-logger option")
elseif(HURRICANE_LOGGER MATCHES "LESS")
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HUR_LESS_LOGGER)
    message(STATUS "The OpenHurricane program is using less-logger option")
elseif(HURRICANE_LOGGER MATCHES "FULL")
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF HUR_FULL_LOGGER)
    message(STATUS "The OpenHurricane program is using full-logger option")
else()
    message(FATAL_ERROR "Invalid logger type")
endif()

#======================================================
# 15. Include source directory
#    Should be add at the end of this file.
#======================================================
add_subdirectory(thirdParty)
add_subdirectory(src)
add_subdirectory(applications)