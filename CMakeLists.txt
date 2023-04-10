cmake_minimum_required(VERSION 3.11 FATAL_ERROR)
project(mashmap VERSION 3.0.0)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Falling back to different standard it not allowed.
set(CMAKE_CXX_EXTENSIONS OFF)  # Make sure no compiler-specific features are used.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)


#set(CMAKE_BUILD_TYPE Release)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: Release Debug Generic." FORCE)
endif()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

option(OPTIMIZE_FOR_NATIVE "Build with -march=native" ON)

if (${CMAKE_BUILD_TYPE} MATCHES Release)
    set(EXTRA_FLAGS "-Ofast -flto -fno-fat-lto-objects")
    if (OPTIMIZE_FOR_NATIVE)
      set(EXTRA_FLAGS "${EXTRA_FLAGS} -march=native")
    endif ()
    set(CMAKE_CXX_FLAGS_RELEASE "-g -DNDEBUG") # reset CXX_FLAGS to replace -O3 with -Ofast
endif ()

if (${CMAKE_BUILD_TYPE} MATCHES Debug)
    # Debug use the defaults - so commenting out:
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O -g")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O -g")
else()
    # Use all standard-compliant optimizations - always add these:
    set (CMAKE_C_FLAGS "${OpenMP_C_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
    set (CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} ${PIC_FLAG} ${EXTRA_FLAGS}")
endif ()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_compile_options(-Wfatal-errors)

include(FetchContent)
FetchContent_Declare(
    assert
    URL https://github.com/jeremy-rifkin/libassert/archive/refs/tags/v1.1.tar.gz
)
FetchContent_MakeAvailable(assert)

add_executable(mashmap
  src/common/utils.cpp
  src/map/mash_map.cpp) 

add_executable(mashmap-align
  src/common/utils.cpp
  src/align/align.cpp) 

target_include_directories(mashmap PRIVATE
  src
  src/common
  )

target_include_directories(mashmap-align PRIVATE
  src
  src/common
  )


target_link_libraries(mashmap
  gsl
  gslcblas
  m
  pthread
  rt
  z
  assert
)

target_link_libraries(mashmap-align
  gsl
  gslcblas
  m
  pthread
  rt
  z
  assert
)

install(TARGETS mashmap DESTINATION bin)
install(TARGETS mashmap-align DESTINATION bin)

# version stuff
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/include)
execute_process(COMMAND bash ${CMAKE_SOURCE_DIR}/scripts/generate_git_version.sh ${CMAKE_SOURCE_DIR}/src)
