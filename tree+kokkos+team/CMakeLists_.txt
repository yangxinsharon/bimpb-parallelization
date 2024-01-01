cmake_minimum_required (VERSION 3.16)
project (bimpbTreeKokkos)

find_package(KokkosKernels REQUIRED)

set(SPACK_CXX $ENV{SPACK_CXX})
if(SPACK_CXX)
  message("found spack compiler ${SPACK_CXX}")
  set(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)  
  set(ENV{CXX} ${SPACK_CXX})
endif()



if(NOT DEFINED KokkosKernels_SOURCE_DIR)
 message( FATAL_ERROR "Please specify KokkosKernels_SOURCE_DIR to include necessary test headers." )
endif()

include_directories(${KokkosKernels_SOURCE_DIR}/test_common)

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")


set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/gmres SRCLIST)
add_library(gmres STATIC ${SRCLIST})

add_executable(bimpb_treekokkos.exe main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_link_libraries(bimpb_treekokkos.exe gmres Kokkos::kokkos Kokkos::kokkoskernels)

