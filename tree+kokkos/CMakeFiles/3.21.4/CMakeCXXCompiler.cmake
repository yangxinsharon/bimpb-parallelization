set(CMAKE_CXX_COMPILER "/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/bin/g++")
set(CMAKE_CXX_COMPILER_ARG1 "")
set(CMAKE_CXX_COMPILER_ID "GNU")
set(CMAKE_CXX_COMPILER_VERSION "9.2.0")
set(CMAKE_CXX_COMPILER_VERSION_INTERNAL "")
set(CMAKE_CXX_COMPILER_WRAPPER "")
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT "14")
set(CMAKE_CXX_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters;cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates;cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates;cxx_std_17;cxx_std_20")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_std_98;cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_std_11;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "cxx_std_14;cxx_aggregate_default_initializers;cxx_attribute_deprecated;cxx_binary_literals;cxx_contextual_conversions;cxx_decltype_auto;cxx_digit_separators;cxx_generic_lambdas;cxx_lambda_init_captures;cxx_relaxed_constexpr;cxx_return_type_deduction;cxx_variable_templates")
set(CMAKE_CXX17_COMPILE_FEATURES "cxx_std_17")
set(CMAKE_CXX20_COMPILE_FEATURES "cxx_std_20")
set(CMAKE_CXX23_COMPILE_FEATURES "")

set(CMAKE_CXX_PLATFORM_ID "Linux")
set(CMAKE_CXX_SIMULATE_ID "")
set(CMAKE_CXX_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_CXX_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_CXX_COMPILER_AR "/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_CXX_COMPILER_RANLIB "/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/bin/gcc-ranlib")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCXX 1)
set(CMAKE_CXX_COMPILER_LOADED 1)
set(CMAKE_CXX_COMPILER_WORKS TRUE)
set(CMAKE_CXX_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_CXX_COMPILER_ENV_VAR "CXX")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_CXX_COMPILER_ID_RUN 1)
set(CMAKE_CXX_SOURCE_FILE_EXTENSIONS C;M;c++;cc;cpp;cxx;m;mm;mpp;CPP;ixx;cppm)
set(CMAKE_CXX_IGNORE_EXTENSIONS inl;h;hpp;HPP;H;o;O;obj;OBJ;def;DEF;rc;RC)

foreach (lang C OBJC OBJCXX)
  if (CMAKE_${lang}_COMPILER_ID_RUN)
    foreach(extension IN LISTS CMAKE_${lang}_SOURCE_FILE_EXTENSIONS)
      list(REMOVE_ITEM CMAKE_CXX_SOURCE_FILE_EXTENSIONS ${extension})
    endforeach()
  endif()
endforeach()

set(CMAKE_CXX_LINKER_PREFERENCE 30)
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES 1)

# Save compiler ABI information.
set(CMAKE_CXX_SIZEOF_DATA_PTR "8")
set(CMAKE_CXX_COMPILER_ABI "ELF")
set(CMAKE_CXX_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_CXX_LIBRARY_ARCHITECTURE "")

if(CMAKE_CXX_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_CXX_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_CXX_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_CXX_COMPILER_ABI}")
endif()

if(CMAKE_CXX_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_CXX_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_CXX_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES "/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/openmpi-3.1.6-liesbt263u4ivo2sdue2ssicegl6zwvo/include;/cm/shared/apps/slurm/current/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/hwloc-1.11.11-vznlzclew7hqyqzj6bsuv5bmsvqq4fva/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/numactl-2.0.14-zk5ignpgj6kb63qlgj2qtbzhoa6uyec7/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libpciaccess-0.16-dv7775isicgotven5m3wodqkamfvz7hu/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cuda-11.1.0-t3z5kvfalh6n7rizd7scrtgskvaz7dxa/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libxml2-2.9.10-robshg56l56shnv2dff3tflbctkhk2a6/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/xz-5.2.5-3doup2tgohwhlvttfm5c2dgyasmtsgti/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libiconv-1.16-ryqlaphegxqe2cjzcw2innllj72ndsti/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/openssl-1.1.1i-3pqsxdfr257u7v5z6vtajzkfp5ngoggp/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/zlib-1.2.11-25wcj7mvn66umzo2nwgpjgpu5l4cxgid/include;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/ncurses-6.2-stjseodsgay3i6wnahbfgzyuqjrljtkw/include;/hpc/applications/nvidia/hpc_sdk/2022_22.2/Linux_x86_64/22.2/compilers/extras/qd/include/qd;/hpc/applications/nvidia/hpc_sdk/2022_22.2/Linux_x86_64/22.2/comm_libs/nvshmem/include;/hpc/applications/nvidia/hpc_sdk/2022_22.2/Linux_x86_64/22.2/comm_libs/nccl/include;/hpc/applications/nvidia/hpc_sdk/2022_22.2/Linux_x86_64/22.2/comm_libs/mpi/include;/hpc/applications/nvidia/hpc_sdk/2022_22.2/Linux_x86_64/22.2/math_libs/include;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/include/c++/9.2.0;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/include/c++/9.2.0/x86_64-pc-linux-gnu;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/include/c++/9.2.0/backward;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/lib/gcc/x86_64-pc-linux-gnu/9.2.0/include;/usr/local/include;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/include;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/lib/gcc/x86_64-pc-linux-gnu/9.2.0/include-fixed;/usr/include")
set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_CXX_IMPLICIT_LINK_DIRECTORIES "/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/kokkos-3.6.01-qu45u5vcjttvscfh76nijf2keukcj5z7/lib64;/cm/shared/apps/slurm/current/lib64;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/lib64;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/lib/gcc/x86_64-pc-linux-gnu/9.2.0;/lib64;/usr/lib64;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/openmpi-3.1.6-liesbt263u4ivo2sdue2ssicegl6zwvo/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/hwloc-1.11.11-vznlzclew7hqyqzj6bsuv5bmsvqq4fva/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/numactl-2.0.14-zk5ignpgj6kb63qlgj2qtbzhoa6uyec7/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libpciaccess-0.16-dv7775isicgotven5m3wodqkamfvz7hu/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/cuda-11.1.0-t3z5kvfalh6n7rizd7scrtgskvaz7dxa/lib64;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libxml2-2.9.10-robshg56l56shnv2dff3tflbctkhk2a6/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/xz-5.2.5-3doup2tgohwhlvttfm5c2dgyasmtsgti/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/libiconv-1.16-ryqlaphegxqe2cjzcw2innllj72ndsti/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/openssl-1.1.1i-3pqsxdfr257u7v5z6vtajzkfp5ngoggp/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/zlib-1.2.11-25wcj7mvn66umzo2nwgpjgpu5l4cxgid/lib;/hpc/spack/opt/spack/linux-centos7-broadwell/gcc-9.2.0/ncurses-6.2-stjseodsgay3i6wnahbfgzyuqjrljtkw/lib;/hpc/spack/opt/spack/linux-centos7-x86_64/gcc-7.3.0/gcc-9.2.0-6zgrndxveon2m5mjhltrqccdcewrdktx/lib")
set(CMAKE_CXX_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
