JULIA_DETAILS=$(julia --version)
CXXWRAP_DETAILS=$(julia --print 'import Pkg; Pkg.status("CxxWrap")')

# --------------- USER INPUTS --------------------------------------------------
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

# --------------- COMPILE CODE -------------------------------------------------
THIS_DIR=$(pwd)
SRC_DIR=deps
COMPILE_DIR=build
SAVE_DIR=src

echo "Removing existing build"
rm -rf $COMPILE_DIR
rm -f $SAVE_DIR/fmm.so

echo "Copying files"
mkdir $COMPILE_DIR
cp -r $SRC_DIR/* $COMPILE_DIR/

echo "Configuring build"
cd $COMPILE_DIR/
./configure
# ./configure --enable-single

echo "Compiling 3d"
cd 3d
# make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB

# -ffast-math might be faster, but not safe in some architectures
make JULIA_H=$JULIA_H JLCXX_H=$JLCXX_H JULIA_LIB=$JULIA_LIB JLCXX_LIB=$JLCXX_LIB EXTRAOBJFLAGS=-ffast-math

cd $THIS_DIR
cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.so

echo -e "\nDone!"

echo -e "\nCompile Summary:"
echo $JULIA_DETAILS
echo $CXXWRAP_DETAILS | sed 's/.*C/C/'

# ---      UNCOMMENT THIS SECTION FOR        ---
# --- ADDITIONAL STEPS FOR BYU SUPERCOMPUTER ---
# rm -f $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.so
#
# cd $COMPILE_DIR/3d/
#
# # libmpi.so often causes trouble and has to be included with fmm
# MPI_LIB=$(dirname $(which mpicxx))/../lib
#
# # Compiler flags with debugging
# # mpicxx -ffast-math -funroll-loops -fabi-version=6 -Wfatal-errors -fopenmp -g -O2 -o fmm fmm-fmm.o -L$JLCXX_LIB -lcxxwrap_julia -fPIC -march=broadwell -Wunused-parameter -Wextra -Wreorder -std=gnu++1z -O3 -DNDEBUG -shared -Wl,-rpath,$MPI_LIB: -LJLCXX_LIB -lcxxwrap_julia -L$JULIA_LIB -ljulia
#
# # Compiler flags without debugging
# mpicxx -ffast-math -funroll-loops -fabi-version=6 -fopenmp -O2 -o fmm fmm-fmm.o -L$JLCXX_LIB -lcxxwrap_julia -fPIC -march=broadwell -std=gnu++1z -O3 -shared -Wl,-rpath,$MPI_LIB: -LJLCXX_LIB -lcxxwrap_julia -L$JULIA_LIB -ljulia
#
# cd $THIS_DIR
# cp $COMPILE_DIR/3d/fmm $SAVE_DIR/fmm.so
#
# # Testing FLOWExaFMM import
# julia -e "import FLOWExaFMM" && echo "FLOWExaFMM installation successful!"
