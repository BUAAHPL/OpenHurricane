# ---------------------------------------------------------------
# Programmer(s): Rao Sihang
# ---------------------------------------------------------------

cmake_minimum_required (VERSION 3.18)
message ( STATUS "The current cmake file is ${CMAKE_CURRENT_LIST_FILE}" ) 

findSourceFiles(${CMAKE_CURRENT_SOURCE_DIR})

if(CUDA_ENABLE)
	findCUDASourceFiles(${CMAKE_CURRENT_SOURCE_DIR})
	get_property(tmp_cuda_src_list GLOBAL PROPERTY global_cuda_src_list) 
	get_property(tmp_cuda_header_list GLOBAL PROPERTY global_cuda_header_list) 
	get_property(tmp_cuda_include_dirs GLOBAL PROPERTY global_cuda_include_dirs) 
	get_property(PRJ_CUDA_tmp_LIBRARIES GLOBAL PROPERTY CUDA_PRJ_LIBRARIES) 
	get_property(CUDA_PRJ_COMPILE_DEF GLOBAL PROPERTY PRJ_COMPILE_DEF) 
	get_property(tmp_header_list GLOBAL PROPERTY global_header_list) 
	get_property(tmp_include_dirs GLOBAL PROPERTY global_include_dirs1) 
	
	list ( APPEND tmp_header_list ${tmp_cuda_header_list} )
	list ( APPEND tmp_include_dirs ${tmp_cuda_include_dirs} )
	list ( APPEND cudaSrcFiles  ${tmp_cuda_src_list} ) 
	
	set(CUDA_LIB "HurricaneCUDALibs")
	add_library( ${CUDA_LIB} STATIC 
		${cudaSrcFiles} ${tmp_header_list}
	)
	target_include_directories ( ${CUDA_LIB}
		PUBLIC 
			${tmp_include_dirs}
	)	
	
	target_compile_definitions ( ${CUDA_LIB}
		PUBLIC
		   ${CUDA_PRJ_COMPILE_DEF} 
	)	
	
	set_target_properties ( ${CUDA_LIB}
		PROPERTIES
			CUDA_SEPARABLE_COMPILATION ON
	)
	set_property(TARGET ${CUDA_LIB} PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
	set_property(GLOBAL APPEND PROPERTY PRJ_LIBRARIES ${CUDA_LIB})
	install(TARGETS ${CUDA_LIB} DESTINATION lib)
endif(CUDA_ENABLE)