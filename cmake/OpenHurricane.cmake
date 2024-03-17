# ---------------------------------------------------------------
# Programmer(s): Rao Sihang
# ---------------------------------------------------------------

function (findSourceFiles currentDir)
	message (STATUS "Finding C/C++ source files from ${currentDir}")
	
	file ( GLOB headerFiles "${currentDir}/*.hpp" "${currentDir}/*.h")
	file ( GLOB sourceFiles "${currentDir}/*.cpp" "${currentDir}/*.c")	
	
	# Get the number of header files
	list ( LENGTH headerFiles nHeaderFiles )
	if ( NOT ( ${nHeaderFiles} EQUAL 0 ) )
		list ( APPEND totalIncludeDirs ${currentDir} )
	endif()
	
	set ( totalHeaderFiles  "" )
	set ( totalSourceFiles  "" )
	set ( totalIncludeDirs  "" )
	
	source_group ( "${currentDir}" FILES ${headerFiles} )
	source_group ( "${currentDir}" FILES ${sourceFiles} )
    list ( APPEND totalHeaderFiles ${headerFiles} )
	list ( APPEND totalSourceFiles  ${sourceFiles} ) 

	file ( GLOB_RECURSE foundFiles LIST_DIRECTORIES true RELATIVE ${currentDir} * )
	foreach (subFound ${foundFiles})
		set (subDir ${currentDir}/${subFound})
		if(IS_DIRECTORY ${subDir})
			file ( GLOB headerFiles "${subDir}/*.hpp" "${subDir}/*.h")
			file ( GLOB sourceFiles "${subDir}/*.cpp" "${subDir}/*.c")
			# Get the number of header files
			list ( LENGTH headerFiles nHeaderFiles )
			if ( NOT ( ${nHeaderFiles} EQUAL 0 ) )
				list ( APPEND totalIncludeDirs ${subDir} )
			endif()
	
			source_group ( "${subDir}" FILES ${headerFiles} )
			source_group ( "${subDir}" FILES ${sourceFiles}    )
			list ( APPEND totalHeaderFiles ${headerFiles} )
			list ( APPEND totalSourceFiles  ${sourceFiles} ) 
		endif()
	endforeach()   
    set_property(GLOBAL APPEND PROPERTY global_src_list "${totalSourceFiles}")
    set_property(GLOBAL APPEND PROPERTY global_header_list "${totalHeaderFiles}")
    set_property(GLOBAL APPEND PROPERTY global_include_dirs1 "${totalIncludeDirs}")
endfunction(findSourceFiles currentDir)

function (findCUDASourceFiles currentDir)
	message (STATUS "Finding CUDA source files from ${currentDir}")

	file ( GLOB headerFiles "${currentDir}/*.hpp" "${currentDir}/*.h" "${currentDir}/*.cuh")
	file ( GLOB sourceFiles "${currentDir}/*.cu" )	
	
	# Get the number of header files
	list ( LENGTH headerFiles nHeaderFiles )
	if ( NOT ( ${nHeaderFiles} EQUAL 0 ) )
		list ( APPEND totalIncludeDirs ${currentDir} )
	endif()

	set ( totalHeaderFiles  "" )
	set ( totalSourceFiles  "" )
	set ( totalIncludeDirs  "" )

	source_group ( "${currentDir}" FILES ${headerFiles} )
	source_group ( "${currentDir}" FILES ${sourceFiles} )
    list ( APPEND totalHeaderFiles ${headerFiles} )
	list ( APPEND totalSourceFiles  ${sourceFiles} ) 

	file ( GLOB_RECURSE foundFiles LIST_DIRECTORIES true RELATIVE ${currentDir} * )
	foreach (subFound ${foundFiles})
		set (subDir ${currentDir}/${subFound})
		if(IS_DIRECTORY ${subDir})
			file ( GLOB headerFiles "${currentDir}/*.hpp" "${currentDir}/*.h" "${currentDir}/*.cuh")
			file ( GLOB sourceFiles "${subDir}/*.cu")
			# Get the number of header files
			list ( LENGTH headerFiles nHeaderFiles )
			if ( NOT ( ${nHeaderFiles} EQUAL 0 ) )
				list ( APPEND totalIncludeDirs ${subDir} )
			endif()

			source_group ( "${subDir}" FILES ${headerFiles} )
			source_group ( "${subDir}" FILES ${sourceFiles}    )
			list ( APPEND totalHeaderFiles ${headerFiles} )
			list ( APPEND totalSourceFiles  ${sourceFiles} ) 
		endif()
	endforeach()   
    set_property(GLOBAL APPEND PROPERTY global_cuda_src_list "${totalSourceFiles}")
    set_property(GLOBAL APPEND PROPERTY global_cuda_header_list "${totalHeaderFiles}")
    set_property(GLOBAL APPEND PROPERTY global_cuda_include_dirs "${totalIncludeDirs}")
endfunction(findCUDASourceFiles currentDir)

macro(hurricane_option OptionName OptionType DocString DEFAULT_VALUE)
   
    if(DEFINED ${OptionName})# if the variable is already defined, keep its value
        set(${OptionName} ${${OptionName}} CACHE ${OptionType} ${DocString} FORCE)       
    else()# if the variable is not already defined, use the default value
        set(${OptionName} ${DEFAULT_VALUE} CACHE ${OptionType} ${DocString})
    endif()

    string(TOUPPER ${${OptionName}} ${OptionName})
endmacro()