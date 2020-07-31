# @Author: Eduardo Alvarez <user>
# @Date:   2017-08-24T12:53:31-06:00
# @Email:  Edo.AlvarezR@gmail.com
# @Last modified by:   user
# @Last modified time: 2018-08-10T20:18:22-06:00

COMPILE_DIR=build
THIS_DIR=$(pwd)
VPM_DIR=/home/edoalvar/Dropbox/FLOWResearch/MyCodes/MyVPM/

echo "Removing existing build"
rm -rf $COMPILE_DIR
rm fmm.so

echo "Copying files"
mkdir $COMPILE_DIR
cp -r src/* $COMPILE_DIR/
# rm -r $COMPILE_DIR/exafmm/3d?*

echo "Configuring build"
cd $COMPILE_DIR/
./configure

echo "Compiling 3d"
cd 3d
make install

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm ./fmm.so

echo "(Replacing fmm.so at $VPM_DIR/src/tools/)"
cp $COMPILE_DIR/3d/fmm $VPM_DIR/src/tools/fmm.so

echo "Done!"
