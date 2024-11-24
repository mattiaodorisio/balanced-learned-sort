cmake_minimum_required(VERSION 3.13)
project(balanced-learned-sort)

# Compile IPS4o disabling parallelism and before setting the flags (it has its own)
option(IPS4O_DISABLE_PARALLEL "Disable parallel code of IPS4o" ON)
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/ips4o)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pthread -O0 -Wall -g -ggdb -fopt-info-vec-optimized")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -pthread -O3 -march=native")
set(CMAKE_BUILD_TYPE RELEASE)

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/extern/learned-sort/third_party/radix_sort)

include_directories(${PROJECT_SOURCE_DIR}/extern/learned-sort/include)
include_directories(${PROJECT_SOURCE_DIR}/extern/pdqsort)
include_directories(${PROJECT_SOURCE_DIR}/extern/ska_sort)

# BLS libraries
include_directories(${PROJECT_SOURCE_DIR}/include)

# Synthetic benchmarks
set(BENCH_SYNTH_uls600 bench_synth_uls600)
add_executable(${BENCH_SYNTH_uls600} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_uls600} PRIVATE ULS600)

set(BENCH_SYNTH_uls1000 bench_synth_uls1000)
add_executable(${BENCH_SYNTH_uls1000} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_uls1000} PRIVATE ULS1000)

set(BENCH_SYNTH_bls bench_synth_bls)
add_executable(${BENCH_SYNTH_bls} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_bls} PRIVATE BLS)

set(BENCH_SYNTH_ls21 bench_synth_ls21)
add_executable(${BENCH_SYNTH_ls21} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_ls21} PRIVATE LS21)

set(BENCH_SYNTH_ips bench_synth_ips)
add_executable(${BENCH_SYNTH_ips} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_ips} PRIVATE IPS)
target_link_libraries(${BENCH_SYNTH_ips} PRIVATE ips4o)

set(BENCH_SYNTH_ls20 bench_synth_ls20)
add_executable(${BENCH_SYNTH_ls20} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_ls20} PRIVATE LS20)

set(BENCH_SYNTH_pdq bench_synth_pdq)
add_executable(${BENCH_SYNTH_pdq} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_pdq} PRIVATE PDQ)

set(BENCH_SYNTH_rs bench_synth_rs)
add_executable(${BENCH_SYNTH_rs} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_rs} PRIVATE RS)
target_link_libraries(${BENCH_SYNTH_rs} PRIVATE radix_sort)

set(BENCH_SYNTH_ska bench_synth_ska)
add_executable(${BENCH_SYNTH_ska} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_ska} PRIVATE SKA)

set(BENCH_SYNTH_stdsort bench_synth_stdsort)
add_executable(${BENCH_SYNTH_stdsort} src/sort_synth.cpp)
target_compile_definitions(${BENCH_SYNTH_stdsort} PRIVATE STDSORT)

################################################################

# Synthetic grid benchmarks
#ULS
set(BENCH_GRID_uls bench_grid_uls)
add_executable(${BENCH_GRID_uls} src/sort_grid.cpp)
target_compile_definitions(${BENCH_GRID_uls} PRIVATE ULS)

#BLS
set(BENCH_GRID_bls bench_grid_bls)
add_executable(${BENCH_GRID_bls} src/sort_grid.cpp)
target_compile_definitions(${BENCH_GRID_bls} PRIVATE BLS)

#LS21
set(BENCH_GRID_ls21 bench_grid_ls21)
add_executable(${BENCH_GRID_ls21} src/sort_grid.cpp)
target_compile_definitions(${BENCH_GRID_ls21} PRIVATE LS21)

################################################################

# Real benchmarks
set(BENCH_REAL_uls600 bench_real_uls600)
add_executable(${BENCH_REAL_uls600} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_uls600} PRIVATE ULS600)

set(BENCH_REAL_uls1000 bench_real_uls1000)
add_executable(${BENCH_REAL_uls1000} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_uls1000} PRIVATE ULS1000)

set(BENCH_REAL_bls bench_real_bls)
add_executable(${BENCH_REAL_bls} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_bls} PRIVATE BLS)

set(BENCH_REAL_ls21 bench_real_ls21)
add_executable(${BENCH_REAL_ls21} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_ls21} PRIVATE LS21)

set(BENCH_REAL_ips bench_real_ips)
add_executable(${BENCH_REAL_ips} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_ips} PRIVATE IPS)
target_link_libraries(${BENCH_REAL_ips} PRIVATE ips4o)

set(BENCH_REAL_ls20 bench_real_ls20)
add_executable(${BENCH_REAL_ls20} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_ls20} PRIVATE LS20)

set(BENCH_REAL_pdq bench_real_pdq)
add_executable(${BENCH_REAL_pdq} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_pdq} PRIVATE PDQ)

set(BENCH_REAL_rs bench_real_rs)
add_executable(${BENCH_REAL_rs} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_rs} PRIVATE RS)
target_link_libraries(${BENCH_REAL_rs} PRIVATE radix_sort)

set(BENCH_REAL_ska bench_real_ska)
add_executable(${BENCH_REAL_ska} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_ska} PRIVATE SKA)

set(BENCH_REAL_stdsort bench_real_stdsort)
add_executable(${BENCH_REAL_stdsort} src/sort_real.cpp)
target_compile_definitions(${BENCH_REAL_stdsort} PRIVATE STDSORT)

# Example
#set(EXAMPLE example)
#add_executable(${EXAMPLE} example/example.cpp)
