cmake_minimum_required (VERSION 3.10)
project (bimpbKokkos)
 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(SPACK_CXX $ENV{SPACK_CXX})
if(SPACK_CXX)
  message("found spack compiler ${SPACK_CXX}")
  set(CMAKE_CXX_COMPILER ${SPACK_CXX} CACHE STRING "the C++ compiler" FORCE)
  set(ENV{CXX} ${SPACK_CXX})
endif()

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/gmres SRCLIST)
add_library(gmres STATIC ${SRCLIST})

add_executable(bimpb_kokkos.exe main_kokkos.cpp readin.cpp matvec_kokkos.cpp pp_timer.cpp)
target_link_libraries(bimpb_kokkos.exe gmres Kokkos::kokkos)

