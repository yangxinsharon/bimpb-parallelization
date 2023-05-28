#ifndef ENV_KOKKOS
#define ENV_KOKKOS

    typedef Kokkos::Serial   HostExecSpace;
    typedef Kokkos::Cuda     DevExecSpace;
    typedef Kokkos::CudaSpace    MemSpace;
    typedef Kokkos::LayoutRight  Layout;
    typedef Kokkos::RangePolicy<HostExecSpace>  host_range_policy;
  	typedef Kokkos::RangePolicy<DevExecSpace>   dev_range_policy;

  	typedef Kokkos::View<double*, Layout, MemSpace>   ViewVectorType;
    typedef Kokkos::View<double**, Layout, MemSpace>  ViewMatrixType;



#endif