export JlCxx_DIR=/Users/randerson/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc/
THIS_DIR=$(pwd)
BUILD_DIR=${THIS_DIR}/build
rm -rf ${BUILD_DIR}
mkdir ${BUILD_DIR}
cd ${BUILD_DIR}
echo "Configuring build..."
echo "  Source path: ${THIS_DIR}/deps/3d"
CC=clang CXX=clang++ cmake ${THIS_DIR}/deps/3d
echo "Making files..."
cmake --build .
echo "Copying dynamic library..."
cp lib/libfmm.dylib ../src/fmm.dylib
