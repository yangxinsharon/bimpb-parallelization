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
typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorInt;
typedef Kokkos::View<double**, Layout, MemSpace>  ViewMatrixInt;
#endif