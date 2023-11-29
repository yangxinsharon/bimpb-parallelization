cmake_minimum_required (VERSION 3.12)
project (bimpbTreeKokkos)

set(SPACK_CXX g++)
if(SPACK_CXX)
  message("found spack compiler g++")
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "the C++ compiler" FORCE)
  set(ENV{CXX} g++)
endif()

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)
find_package(KokkosKernels REQUIRED)

if(NOT DEFINED KokkosKernels_SOURCE_DIR)
 message( FATAL_ERROR "Please specify KokkosKernels_SOURCE_DIR to include necessary test headers." )
endif()

include_directories(${KokkosKernels_SOURCE_DIR}/test_common)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/gmres SRCLIST)
add_library(gmres STATIC ${SRCLIST})

add_executable(bimpb_treekokkos.exe main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_link_libraries(bimpb_treekokkos.exe gmres Kokkos::kokkos Kokkos::kokkoskernels)
