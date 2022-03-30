# --------------- SUPERCOMPUTER ------------------------------------------------
# module load julia/1.6
module load gcc/10
module load openmpi/4.1

# --------------- USER INPUTS --------------------------------------------------
# choose your C compiler (use LLVM on MacOS)
# CC=/apps/gcc/10.2.0/bin/gcc
CC=/apps/openmpi/4.1.1/gcc-10.2.0_cuda-11.2.1/bin/mpicc

# choose your C++ compiler (use LLVM on MacOS)
# CXX=/apps/gcc/10.2.0/bin/g++
#CXX=/bin/g++
CXX=/apps/openmpi/4.1.1/gcc-10.2.0_cuda-11.2.1/bin/mpixx

# JULIA_H must point to the directory that contains julia.h
# JULIA_H=/fslhome/rander39/julia/julia-1.7.2/include/julia/
JULIA_H=/apps/julia/1.6.1/gcc-10.2.0/include/julia

# JLCXX_H must point to the directory that contains jlcxx/jlcxx.hpp from CxxWrap
# NOTE: You can find this by typing `CxxWrap.prefix_path()` in the Julia REPL
# JLCXX_H=/fslhome/rander39/.julia/artifacts/6fbda767b0cb63d5cd1b258fd698b8fb14d655c0/include/
JLCXX_H=/fslhome/rander39/.julia/artifacts/4fcd159fccd2f12b8c8c3da884709dc1de7a30ae/include/

# Julia_LIB must point to the directory that contains libjulia.so.x
JULIA_LIB=$JULIA_H/../../lib

# JLCXX_LIB must point to the directory that contains libcxxwrap_julia.so.0.x.x
JLCXX_LIB=$JLCXX_H/../lib

# MPICXX path
# MPIHOME=/usr/local/bin/mpicxx

# OpenMP flags
# LDFLAGS="-L/usr/local/opt/llvm/lib"
# -Wl,-rpath,/usr/local/opt/llvm/lib"

# NOTE: on mac, shared libraries have the .dylib extension. This may require further configuration depending on how Eduardo has set this up.

THIS_DIR="$HOME/julia/dev/FLOWExaFMM"

# --------------- COMPILE CODE -------------------------------------------------
SRC_DIR=$THIS_DIR/deps
COMPILE_DIR=$THIS_DIR/build
SAVE_DIR=$THIS_DIR/src/fmm_tmp

echo "Copying files"
mkdir $COMPILE_DIR
cp -r $SRC_DIR/* $COMPILE_DIR/

echo "Configuring build"
cd $COMPILE_DIR/
sh $COMPILE_DIR/configure CXX=$CXX CC=$CC --prefix=$COMPILE_DIR/3d
#HOME/.julia/dev/FLOWExaFMM/deps/3d
# ./configure CXX=$CXX MPICXX=$MPIHOME LDFLAGS=$LDFLAGS --enable-single
# CC=$CC
# ./configure

echo "Compiling 3d"
cd 3d
make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB

tmp_name=$(mktemp $SAVE_DIR/fmm.XXXXXXXX)
cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $tmp_name
export FMM=$tmp_name

echo "Done!"
