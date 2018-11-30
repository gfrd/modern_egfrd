
function(check_make_unique)
  
  if (NOT DEFINED HAS_MAKE_UNIQUE)
 
    set (_BINDIR "${CMAKE_BINARY_DIR}/test-make_unique")
    set (_SRCFILE ${CMAKE_CURRENT_LIST_DIR}/cmake/test-make_unique.cpp )
    message (STATUS "Checking support for std::make_unique")
	
    if (NOT EXISTS "${_SRCFILE}")
      message (STATUS "Checking support for std::make_unique: not supported (test not found)")
      set (HAS_MAKE_UNIQUE FALSE CACHE INTERNAL "support for std::make_unique")
      return()
    endif ()

	IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 3.8)
		try_compile (HAS_MAKE_UNIQUE "${_BINDIR}" "${_SRCFILE}")
	ELSE()
		try_compile (HAS_MAKE_UNIQUE "${_BINDIR}" "${_SRCFILE}" CXX_STANDARD 11)
	ENDIF()

    if (HAS_MAKE_UNIQUE)
      message (STATUS "Checking support for std::make_unique - done")
    else()
      message (STATUS "Checking support for std::make_unique - not supported")
    endif ()

    set (HAS_MAKE_UNIQUE ${HAS_MAKE_UNIQUE} CACHE INTERNAL "support for std::make_unique")
  endif ()

endfunction()