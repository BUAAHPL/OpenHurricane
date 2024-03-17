# ---------------------------------------------------------------
# Programmer(s): Rao Sihang
# ---------------------------------------------------------------
# ---------------------------------------------------------------
# Setup the CUDA languge and CUDA libraries.
# ---------------------------------------------------------------

if(NOT CMAKE_CUDA_HOST_COMPILER)
  # If a user did not provide the host compiler, then CXX compiler was used.
  set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "NVCC host compiler")
  message(STATUS "Setting Cmake CUDA host compiler: ${CMAKE_CUDA_HOST_COMPILER}")
endif()

if(${CMAKE_VERSION} VERSION_LESS "3.18.0")
    set(CMAKE_CUDA_STANDARD 14)
else()
    set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})
endif()

if(POLICY CMP0146)
  cmake_policy(SET CMP0146 OLD)
endif()

if(POLICY CMP0148)
  cmake_policy(SET CMP0148 OLD)
endif()
find_package(CUDA REQUIRED)
if(CUDA_FOUND)
	set ( CMAKE_CUDA_STANDARD_REQUIRED ON )
	enable_language(CUDA)
    set_property(GLOBAL APPEND PROPERTY PRJ_COMPILE_DEF CUDA_PARALLEL)
else(CUDA_FOUND)        
    message(FATAL_ERROR "The CUDA was not found.")
endif(CUDA_FOUND)



if(${CUDA_VERSION_STRING} VERSION_LESS "12.1")
	set ( CMAKE_CUDA_ARCHITECTURES "35;50;52;60;61;70;72;75;80;86;87" )
else()
	set ( CMAKE_CUDA_ARCHITECTURES "50;52;60;61;70;72;75;80;86;87" )
endif()

# Set CUDA standard
set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD})

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

# Make the CUDA_rt_LIBRARY advanced like the other CUDA_*_LIBRARY variables
mark_as_advanced(FORCE CUDA_rt_LIBRARY)

# Show CUDA flags
mark_as_advanced(CLEAR CMAKE_CUDA_FLAGS)

# Set CUDA flags
set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Wno-deprecated-gpu-targets" )
set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda" )
set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-extended-lambda")
if(CMAKE_BUILD_TYPE STREQUAL "Debug" )
   set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -g -m64" )
else(CMAKE_BUILD_TYPE STREQUAL "Debug" ) # Release
   set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -O3 -m64" )
endif(CMAKE_BUILD_TYPE STREQUAL "Debug" )
 
set(CUDA_LIBRARIES "${CUDA_LIBRARIES}" "${CUDA_cublas_LIBRARY}")
set(CUDA_LIBRARIES "${CUDA_LIBRARIES}" "${CUDA_cusparse_LIBRARY}")
set(CUDA_LIBRARIES "${CUDA_LIBRARIES}" "${CUDA_cusolver_LIBRARY}")
set(CUDA_LIBRARIES "${CUDA_LIBRARIES}" "${CUDA_cudadevrt_LIBRARY}")
 
set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES "${CUDA_LIBRARIES}")
set_property(GLOBAL APPEND PROPERTY CUDA_PRJ_LIBRARIES "${CUDA_LIBRARIES}")

# Print out information about CUDA.
message(STATUS "The CUDA will be used")
message(STATUS "CUDA Version:               ${CUDA_VERSION_STRING}")
message(STATUS "CUDA Architectures:         ${CMAKE_CUDA_ARCHITECTURES}")
message(STATUS "CUDA Compiler:              ${CMAKE_CUDA_COMPILER}")
message(STATUS "CUDA Host Compiler:         ${CMAKE_CUDA_HOST_COMPILER}")
message(STATUS "CUDA Include Path:          ${CUDA_INCLUDE_DIRS}")
message(STATUS "CUDA Libraries:             ${CUDA_LIBRARIES}")
message(STATUS "CUDA Compile Flags:         ${CMAKE_CUDA_FLAGS}")
message(STATUS "CUDA Link Flags:            ${CMAKE_CUDA_LINK_FLAGS}")
message(STATUS "CUDA Link Executable:       ${CMAKE_CUDA_LINK_EXECUTABLE}")
message(STATUS "CUDA Separable Compilation: ${CMAKE_CUDA_SEPARABLE_COMPILATION}")
message(STATUS "CUDA standard:              ${CMAKE_CUDA_STANDARD}")