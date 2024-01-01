#ifndef GL_CONSTANTS
#define GL_CONSTANTS
const double pi=3.14159265358979324;
const double one_over_4pi=0.079577471545948;
const double bulk_coef=8.430325455;
const double units_coef=332.0716;
const double epsw=80.0; 
const double epsp=1.0;
const double eps=80.0;
const double bulk_strength=0.15;   	//ion_strength in M
const double kappa2=0.0158068602;	//kappa2=bulk_coef*bulk_strength/epsw;
const double kappa=0.1257253365;

const int order=3;
const int maxparnode=100;
const double theta=0.8;
#endif

// add_executable(bimpb_treekokkos.openmp main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
// target_compile_definitions(bimpb_treekokkos.openmp PRIVATE -DUSE_OPENMP)
// target_link_libraries(bimpb_treekokkos.openmp gmres Kokkos::kokkos)

// add_executable(bimpb_treekokkos.cuda main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
// target_compile_definitions(bimpb_treekokkos.cuda PRIVATE -DUSE_CUDA)
// target_link_libraries(bimpb_treekokkos.cuda gmres Kokkos::kokkos)

// add_executable(bimpb_treekokkos.uvm main_kokkos.cpp readin.cpp matvec_kokkos.cpp treecode.cpp utilities.cpp pp_timer.cpp)
// target_compile_definitions(bimpb_treekokkos.uvm PRIVATE -DUSE_UVM)
// target_link_libraries(bimpb_treekokkos.uvm gmres Kokkos::kokkos)