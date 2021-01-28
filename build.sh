# set environment variables to use LLVM compilers and libraries
export CMAKE_C_COMPILER_JLENV=/usr/local/Cellar/llvm/10.0.1/bin/clang
export CMAKE_CXX_COMPILER_JLENV=/usr/local/Cellar/llvm/10.0.1/bin/clang++
export OPENMP_LIBRARIES_JLENV=/usr/local/Cellar/llvm/10.0.1/lib
export OPENMP_INCLUDES_JLENV=/usr/local/Cellar/llvm/10.0.1/include
# set environment variables for use with CxxWrap
export JlCxx_DIR=/Users/randerson/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc

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
make

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.so

echo "Done!"
