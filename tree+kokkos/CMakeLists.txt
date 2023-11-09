cmake_minimum_required (VERSION 3.10)
project (bimpbTreeKokkos)

set(SPACK_CXX g++)
if(SPACK_CXX)
  message("found spack compiler g++")
  set(CMAKE_CXX_COMPILER g++ CACHE STRING "the C++ compiler" FORCE)
  set(ENV{CXX} g++)
endif()

set(Kokkos_DIR "$ENV{Kokkos_ROOT}" CACHE STRING "Kokkos root directory")
find_package(Kokkos REQUIRED)

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/gmres SRCLIST)
add_library(gmres STATIC ${SRCLIST})

add_executable(bimpb_treekokkos.exe main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_link_libraries(bimpb_treekokkos.exe gmres Kokkos::kokkos)

add_executable(bimpb_treekokkos.openmp main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_compile_definitions(bimpb_treekokkos.openmp PRIVATE -DUSE_OPENMP)
target_link_libraries(bimpb_treekokkos.openmp gmres Kokkos::kokkos)

add_executable(bimpb_treekokkos.cuda main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_compile_definitions(bimpb_treekokkos.cuda PRIVATE -DUSE_CUDA)
target_link_libraries(bimpb_treekokkos.cuda gmres Kokkos::kokkos)

add_executable(bimpb_treekokkos.uvm main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
target_compile_definitions(bimpb_treekokkos.uvm PRIVATE -DUSE_UVM)
target_link_libraries(bimpb_treekokkos.uvm gmres Kokkos::kokkos)

