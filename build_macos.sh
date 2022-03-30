# --------------- USER INPUTS --------------------------------------------------
# choose your C compiler (use LLVM on MacOS)
CC=/usr/local/opt/llvm/bin/clang

# choose your C++ compiler (use LLVM on MacOS)
CXX=/usr/local/opt/llvm/bin/clang++

# JULIA_H must point to the directory that contains julia.h
JULIA_H=/Applications/Julia-1.6.app/Contents/Resources/julia/include/julia

# JLCXX_H must point to the directory that contains jlcxx/jlcxx.hpp from CxxWrap
# NOTE: You can find this by typing `CxxWrap.prefix_path()` in the Julia REPL
# JLCXX_H=/Users/randerson/.julia/artifacts/eb1ece2a20e3bd68968d861ee1e8e6a00d077d9e/include
JLCXX_H=/Users/randerson/tmps/juliadev/libcxxwrap-julia/include

# Julia_LIB must point to the directory that contains libjulia.so.x
JULIA_LIB=$JULIA_H/../../lib

# JLCXX_LIB must point to the directory that contains libcxxwrap_julia.so.0.x.x
# JLCXX_LIB=$JLCXX_H/../lib
JLCXX_LIB=/Users/randerson/tmps/juliadev/libcxxwrap-julia/build_llvm_1_6/lib

# MPICXX path
MPIHOME=/usr/local/bin/mpicxx

# OpenMP flags
LDFLAGS="-L/usr/local/opt/llvm/lib"
# -Wl,-rpath,/usr/local/opt/llvm/lib"

# NOTE: on mac, shared libraries have the .dylib extension. This may require further configuration depending on how Eduardo has set this up.

# --------------- COMPILE CODE -------------------------------------------------
THIS_DIR=$(pwd)
SRC_DIR=deps
COMPILE_DIR=build
SAVE_DIR=src

echo "Removing existing build"
rm -rf $COMPILE_DIR
rm $SAVE_DIR/fmm.so

echo "Copying files"
mkdir $COMPILE_DIR
cp -r $SRC_DIR/* $COMPILE_DIR/

echo "Configuring build"
cd $COMPILE_DIR/
./configure CC=$CC CXX=$CXX MPICXX=$MPIHOME LDFLAGS=$LDFLAGS --disable-mpi --disable-omp
# --enable-single
# ./configure

echo "Compiling 3d"
cd 3d
make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.dylib

echo "Done!"