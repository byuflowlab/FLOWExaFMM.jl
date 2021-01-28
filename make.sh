# @

# THIS_DIRECTORY=$(pwd)
SRC_DIR="$THIS_DIRECTORY/deps"
COMPILE_DIR="$THIS_DIRECTORY/build"
SAVE_DIR="$THIS_DIRECTORY/src"

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
