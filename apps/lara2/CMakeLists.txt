# ===========================================================================
#                  SeqAn - The Library for Sequence Analysis
# ===========================================================================
# File: /apps/lara2/CMakeLists.txt
#
# CMakeLists.txt file for lara2
#
# ==========================================================================
# Author: Gianvito Urgese <gianvito.urgese@polito.it>
# ==========================================================================

cmake_minimum_required(VERSION 3.5)
project (seqan_apps_lara2 CXX)
message (STATUS "Configuring apps/lara")

set(CMAKE_CXX_STANDARD  14)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
FIND_PACKAGE(OpenMP REQUIRED)
if(OPENMP_FOUND)
    MESSAGE("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    MESSAGE("OPENMP NOT FOUND")
endif()
#find_package (SeqAn REQUIRED)

# Add include directories.
#
SET(SEQAN_INCLUDE_DIRS "../../include")
include_directories(${SEQAN_INCLUDE_DIRS} ) # FIXME make this variable
include_directories(./include)

set(SOURCE_FILES
        src/laragu.cpp)

add_executable(laragu ${SOURCE_FILES})

# Vienna RNA static library should be installed
find_library(RNA libRNA.a PATHS)
if(RNA)
    message("Vienna RNA FOUND")
else()
    message("Vienna RNA NOT FOUND")
endif()

find_package(Boost COMPONENTS system filesystem REQUIRED)

target_link_libraries(laragu
        ${RNA}
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        )
