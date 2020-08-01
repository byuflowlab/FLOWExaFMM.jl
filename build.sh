COMPILE_DIR=build
THIS_DIR=$(pwd)

echo "Removing existing build"
rm -rf $COMPILE_DIR
rm fmm.so

echo "Copying files"
mkdir $COMPILE_DIR
cp -r src/* $COMPILE_DIR/

echo "Configuring build"
cd $COMPILE_DIR/
./configure

echo "Compiling 3d"
cd 3d
make

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm ./fmm.so

echo "Done!"
