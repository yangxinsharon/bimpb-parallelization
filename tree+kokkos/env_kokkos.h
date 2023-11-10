#ifndef ENV_KOKKOS
#define ENV_KOKKOS

#include <Kokkos_Core.hpp>
typedef Kokkos::Serial   HostExecSpace;
typedef Kokkos::Cuda     DevExecSpace;
typedef Kokkos::CudaSpace    MemSpace;
typedef Kokkos::LayoutRight  Layout;
typedef Kokkos::RangePolicy<HostExecSpace>  host_range_policy;
typedef Kokkos::RangePolicy<DevExecSpace>   dev_range_policy;

typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorDouble;
typedef Kokkos::View<double**, Layout, MemSpace>  ViewMatrixDouble;
typedef Kokkos::View<int*, Layout, MemSpace>   ViewVectorInt;
typedef Kokkos::View<int**, Layout, MemSpace>  ViewMatrixInt;
typedef Kokkos::View<double***, Layout, MemSpace>  TensorDouble;

typedef Kokkos::TeamPolicy<DevExecSpace>  team_policy;
typedef Kokkos::TeamPolicy<DevExecSpace>::member_type  member_type;
#endif

// #if defined(USE_UVM)

//   // set the Kokkos execution space and range policy
// #if defined(USE_OPENMP)
//   typedef Kokkos::OpenMP     ExecSpace;
//   typedef Kokkos::HostSpace  MemSpace;
//   std::cout << "using Kokkos with OpenMP backend:\n";
// #elif defined(USE_CUDA)
//   typedef Kokkos::Cuda       ExecSpace;
//   typedef Kokkos::CudaSpace  MemSpace;
//   std::cout << "using Kokkos with CUDA backend:\n";
// #elif defined(USE_UVM)
//   typedef Kokkos::Cuda          ExecSpace;
//   typedef Kokkos::CudaUVMSpace  MemSpace;
//   std::cout << "using Kokkos with CUDA backend:\n";
// #else
//   typedef Kokkos::Serial     ExecSpace;
//   typedef Kokkos::HostSpace  MemSpace;
//   std::cout << "using Kokkos with Serial backend:\n";
// #endif
//   typedef Kokkos::RangePolicy<ExecSpace>            RangePol;
//   typedef Kokkos::View<double*, MemSpace>           VecViewDev;
//   typedef Kokkos::View<double*, Kokkos::HostSpace>  VecViewHost;
