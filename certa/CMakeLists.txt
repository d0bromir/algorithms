cmake_minimum_required(VERSION 3.18)
project(certa LANGUAGES CXX VERSION 0.1)

# CERTA_CUDA: AUTO (build GPU back-end if nvcc is found), ON, or OFF.
set(CERTA_CUDA AUTO CACHE STRING "Build the CUDA back-end (AUTO/ON/OFF)")
# A100 = 80, RTX 30xx = 86, H100 = 90.
set(CERTA_CUDA_ARCHS "80;86" CACHE STRING "CUDA architectures")
option(CERTA_NATIVE "Tune for the build machine's CPU" ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(CheckCXXCompilerFlag)
if(CERTA_NATIVE)
  # x86-64 and recent aarch64 GCC/Clang accept -march=native; older aarch64
  # compilers only accept -mcpu=native.
  check_cxx_compiler_flag(-march=native HAVE_MARCH_NATIVE)
  check_cxx_compiler_flag(-mcpu=native HAVE_MCPU_NATIVE)
  if(HAVE_MARCH_NATIVE)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-march=native>)
  elseif(HAVE_MCPU_NATIVE)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mcpu=native>)
  endif()
endif()

find_package(Threads REQUIRED)
find_package(ZLIB)

set(USE_CUDA OFF)
if(NOT CERTA_CUDA STREQUAL "OFF")
  include(CheckLanguage)
  check_language(CUDA)
  if(CMAKE_CUDA_COMPILER)
    set(CMAKE_CUDA_ARCHITECTURES ${CERTA_CUDA_ARCHS})
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(USE_CUDA ON)
  elseif(CERTA_CUDA STREQUAL "ON")
    message(FATAL_ERROR "CERTA_CUDA=ON but no CUDA compiler was found")
  endif()
endif()

add_library(certa_core STATIC src/index.cpp src/seqio.cpp src/mapper_cpu.cpp)
target_include_directories(certa_core PUBLIC include src)
target_link_libraries(certa_core PUBLIC Threads::Threads)
if(ZLIB_FOUND)
  target_compile_definitions(certa_core PUBLIC CERTA_WITH_ZLIB)
  target_link_libraries(certa_core PUBLIC ZLIB::ZLIB)
endif()

if(USE_CUDA)
  add_executable(certa src/main.cpp src/mapper_gpu.cu)
  set_target_properties(certa PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)
  target_compile_options(certa PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-O3 --expt-relaxed-constexpr>)
else()
  add_executable(certa src/main.cpp src/gpu_stub.cpp)
endif()
target_link_libraries(certa PRIVATE certa_core)

enable_testing()
add_executable(test_certa tests/test_certa.cpp src/gpu_stub.cpp)
target_link_libraries(test_certa PRIVATE certa_core)
add_test(NAME certificate COMMAND test_certa)

message(STATUS "CERTA: CUDA back-end ${USE_CUDA} (archs ${CERTA_CUDA_ARCHS}), zlib ${ZLIB_FOUND}")
