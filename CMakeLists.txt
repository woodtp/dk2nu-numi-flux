cmake_minimum_required(VERSION 3.21)
project(flux1)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(ROOT REQUIRED)
# include_directories($ENV{DK2NU_INC}/dk2nu/tree ${ROOT_INCLUDE_DIRS})
include_directories($ENV{DK2NU_INC}/dk2nu/tree ${ROOT_INCLUDE_DIRS})
link_directories($ENV{DK2NU_LIB} $ENV{ROOTSYS}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Werror -Wextra")

add_executable(makeFluxHistograms makeFluxHistograms.cc)

target_link_libraries(makeFluxHistograms dk2nuTree ${ROOT_LIBRARIES})
