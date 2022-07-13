# --------------- USER INPUTS --------------------------------------------------
# choose your C++ compiler
CXX=/opt/homebrew/bin/gcc-11

# JULIA_H must point to the directory that contains julia.h
# NOTE: You can find this by typing `abspath(Sys.BINDIR, Base.INCLUDEDIR)` in
#       the Julia REPL
JULIA_H=$(julia --print "abspath(Sys.BINDIR, Base.INCLUDEDIR)")
JULIA_H=${JULIA_H%\"}; JULIA_H=${JULIA_H#\"}; JULIA_H=$JULIA_H/julia

# JLCXX_H must point to the directory that contains jlcxx/jlcxx.hpp from CxxWrap
# NOTE: You can find this by typing `CxxWrap.prefix_path()` in the Julia REPL
JLCXX_H=$(julia --print "import CxxWrap; CxxWrap.prefix_path()")
JLCXX_H=${JLCXX_H%\"}; JLCXX_H=${JLCXX_H#\"}; JLCXX_H=$JLCXX_H/include

# Julia_LIB must point to the directory that contains libjulia.so.x
JULIA_LIB=$JULIA_H/../../lib

# JLCXX_LIB must point to the directory that contains libcxxwrap_julia.so.0.x.x
JLCXX_LIB=$JLCXX_H/../lib

# MPICXX path
MPIHOME=/usr/local/bin/mpicxx

# OpenMP flags
LDFLAGS="-L/opt/homebrew/opt/llvm/lib"
# LDFLAGS="-L/usr/local/opt/llvm/lib"
# -Wl,-rpath,/usr/local/opt/llvm/lib"

# NOTE: on mac, shared libraries have the .dylib extension.

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
./configure CXX=$CXX
# ./configure CXX=$CXX LDFLAGS=$LDFLAGS --disable-mpi
# ./configure CC=$CC CXX=$CXX MPICXX=$MPIHOME LDFLAGS=$LDFLAGS
# ./configure CC=$CC CXX=$CXX MPICXX=$MPIHOME LDFLAGS=$LDFLAGS --disable-mpi --disable-omp
# --enable-single
# ./configure

echo "Compiling 3d"
cd 3d
make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.dylib

echo "Done!"
