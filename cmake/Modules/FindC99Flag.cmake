INCLUDE (${CMAKE_ROOT}/Modules/FindPackageHandleStandardArgs.cmake)
INCLUDE (${CMAKE_ROOT}/Modules/CheckCSourceCompiles.cmake)

# If the compiler flags have already been set, return now
IF(C99_C_FLAGS)
    RETURN()
ENDIF(C99_C_FLAGS)

# Define the code used to test C99
SET(TEST_CODE
"
int main () {
    for (int i = 0; i < 6; i++) 
        ;
}
")

# Define the possible flags
SET(C99_C_FLAG_CANDIDATES
    # Clang and Portland Group (C99 by default)
    " "
    # GNU, Intel
    "-std=c99"
    # Intel Windows
    "/Qstd=c99"
)

# Run over flag options and see what is necessary
SET(C99_FOUND 0)
FOREACH(FLAG ${C99_C_FLAG_CANDIDATES})
    SET(SAFE_CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
    SET(CMAKE_REQUIRED_FLAGS "${FLAG}")
    UNSET(COMPILER_SUPPORTS_C99 CACHE)
    MESSAGE (STATUS "Try C99 C flag = [${FLAG}]")
    CHECK_C_SOURCE_COMPILES("${TEST_CODE}" COMPILER_SUPPORTS_C99)
    SET(CMAKE_REQUIRED_FLAGS "${SAFE_CMAKE_REQUIRED_FLAGS}")
    IF(COMPILER_SUPPORTS_C99)
        SET(C99_FOUND 1)
        SET(C99_FLAG_INTERNAL "${FLAG}")
        BREAK()
    ENDIF(COMPILER_SUPPORTS_C99)
ENDFOREACH(FLAG ${C99_C_FLAG_CANDIDATES})

# Store the variable
IF(C99_FOUND)
    SET(C99_C_FLAGS "${C99_FLAG_INTERNAL}", CACHE
        STRING "C Compiler flag to indicate C99 standard")
ELSE()
    MESSAGE (FATAL_ERROR "This C compiler does not support C99.")
ENDIF(C99_FOUND)

FIND_PACKAGE_HANDLE_STANDARD_ARGS (C99_FLAG DEFAULT_MSG C99_C_FLAGS)

MARK_AS_ADVANCED(C99_C_FLAGS)
