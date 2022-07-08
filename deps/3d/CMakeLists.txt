set(CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER_JLENV}) # "/usr/local/Cellar/llvm/10.0.1/bin/clang")
set(CMAKE_CXX_COMPILER $ENV{CMAKE_CXX_COMPILER_JLENV}) # "/usr/local/Cellar/llvm/10.0.1/bin/clang++")
set(OPENMP_LIBRARIES $ENV{OPENMP_LIBRARIES_JLENV}) # "/usr/local/Cellar/llvm/10.0.1/lib")
set(OPENMP_INCLUDES $ENV{OPENMP_INCLUDES_JLENV}) # "/usr/local/Cellar/llvm/10.0.1/include")

project(FMM CXX)
cmake_minimum_required(VERSION 3.12)  # min required version for os x with openmp
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_VERBOSE_MAKEFILE ON)  # turn on if you want verbose
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
# set(CMAKE_OSX_SYSROOT /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk)
# apple-specific changes
message(STATUS "CMAKE_CXX_COMPILER set as ${CMAKE_CXX_COMPILER}")
set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -Wno-unused-command-line-argument")
set(OpenMP_CXX_LIB_NAMES "omp")
set(OpenMP_omp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
# set traverse flag
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS} -DEXAFMM_EAGER")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DEXAFMM_EAGER")
# auto find cxx and openmp
find_package(JlCxx)
find_package(OpenMP REQUIRED)
get_target_property(JlCxx_location JlCxx::cxxwrap_julia LOCATION)
get_filename_component(JlCxx_location ${JlCxx_location} DIRECTORY)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${JlCxx_location}")
message(STATUS "Found JlCxx at ${JlCxx_location}")
# set include and library locations
add_library(fmm SHARED fmm.cxx)
target_include_directories(fmm PRIVATE ${CMAKE_SOURCE_DIR})
# which linker flags are used?
message(STATUS "CMAKE_EXE_LINKER_FLAGS_INIT set as ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "CMAKE_SHARED_LINKER_FLAGS_INIT set as ${CMAKE_SHARED_LINKER_FLAGS}")
message(STATUS "CMAKE_MODULE_LINKER_FLAGS_INIT set as ${CMAKE_MODULE_LINKER_FLAGS}")
message(STATUS "CMAKE_OPENMP_LINKER_FLAGS_INIT set as ${CMAKE_OPENMP_LINKER_FLAGS_INIT}")
target_link_directories(fmm PRIVATE "/usr/local/lib")  # USER: probably not needed on linux, but not auto found in mac/homebrew
target_link_libraries(fmm PRIVATE JlCxx::cxxwrap_julia OpenMP::OpenMP_CXX)
install(TARGETS
  fmm
LIBRARY DESTINATION lib
ARCHIVE DESTINATION lib
RUNTIME DESTINATION lib)
