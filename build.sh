# --------------- USER INPUTS --------------------------------------------------
# JULIA_H must point to the directory that contains julia.h
JULIA_H=/home/edoalvar/Programs/julia-1.6.1/include/julia

# JLCXX_H must point to the directory that contains jlcxx/jlcxx.hpp from CxxWrap
# NOTE: You can find this by typing `CxxWrap.prefix_path()` in the Julia REPL
JLCXX_H=/home/edoalvar/.julia/artifacts/16e1de4679fb8520a8af4e6831c7c8e9893d18b4/include

# Julia_LIB must point to the directory that contains libjulia.so.x
JULIA_LIB=$JULIA_H/../../lib

# JLCXX_LIB must point to the directory that contains libcxxwrap_julia.so.0.x.x
JLCXX_LIB=$JLCXX_H/../lib


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
# ./configure --enable-single
./configure

echo "Compiling 3d"
cd 3d
make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.so

echo "Done!"
