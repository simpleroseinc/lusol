cmake_minimum_required(VERSION 3.16)

project(lusol LANGUAGES C)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

set(CMAKE_VS_JUST_MY_CODE_DEBUGGING ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# set the sources for the main
file(GLOB_RECURSE lusol_srcs "commonlib.c" "hbio.c" "lusol.c" "lusolio.c" "lusolmain.c" "mmio.c" "myblas.c")

# set the target executable
add_library(lusol SHARED ${lusol_srcs})

# set properties of lusol
set_target_properties(lusol PROPERTIES ARCHIVE_OUTPUT_DIRECTORY lib)

# set the include directories lusol
target_include_directories(lusol PUBLIC
	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)