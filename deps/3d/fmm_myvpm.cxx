/**
 * @Author: Eduardo Alvarez <user>
 * @Date:   2017-08-24T16:10:58-06:00
 * @Email:  Edo.AlvarezR@gmail.com
 * @Last modified by:   user
 * @Last modified time: 2018-09-27T16:38:23-06:00
 * @Comments: I replaced fmm.cxx from the original build to make this module
 *            an interface between MyVPM and exafmm. The original module can be
 *            found at fmm-org.cxx
 */
 // ExaFMM modules
#include "args.h"
#include "build_tree.h"
#include "dataset.h"
#include "kernel.h"
#include "timer.h"
#if EXAFMM_EAGER
#include "traverse_eager.h"
#elif EXAFMM_LAZY
#include "traverse_lazy.h"
#endif
#include "verify.h"
using namespace exafmm;

// Tools for interfacing
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>  //for 'typeid' to work
#include <string>
#include "jlcxx/jlcxx.hpp"  //C++ wrapper for julia

#define M_PIl          3.141592653589793238462643383279502884L

// ---------- MyVPM data structure ---------------------------------------------
const int _n_param = 8; // Number of parameters that define a particle
// Index of each parameter (from 1 to _n_param)
const int _x_i = 1;
const int _y_i = 2;
const int _z_i = 3;
const int _gamma_x_i = 4;
const int _gamma_y_i = 5;
const int _gamma_z_i = 6;
const int _blob_sigma_i = 7;
const int _vol_i = 8;
// States the existence of the parameter
const int _params_i[_n_param] = {_x_i, _y_i, _z_i,
                   _gamma_x_i, _gamma_y_i, _gamma_z_i,
                   _blob_sigma_i, _vol_i};
// -----------------------------------------------------------------------------

const char read_file[32] = "verify_kernel.pf";    // Name of input .pf file for main
const char delimiter = '\t';                      // Delimiter in file
const int header_length = 6;                      // Number of lines in header
int np;                                           // Number of particles
int nparam;                                       // Parameters per particle
std::vector<real_t> probe = {1.0, 1.0, 1.41421};  // Where to probe the potential

#include "fmm_utils.h"

// int main(int argc, char ** argv){
int main(){
  char aux0 = 'a';
  char * aux1 = &aux0;
  int argc = 1;
  char ** argv = &aux1;
  // Reads file of particles
  std::cout << "Reading file\n";
  Bodies particles = readBodies(read_file, true);

  // Adds the point where to probe the potential as a particle
  std::cout << "Adding probe\n";
  Bodies bodies(np+1);
  for(int i=0; i<np; ++i) bodies[i] = particles[i];
  for(int i=0; i<3; ++i) bodies[np].X[i] = probe[i];
  bodies[np].q = 0;

  // Overwrites initialization values of ExaFMM
  std::cout << "Overwriting FMM parameters\n";
  std::cout << "-------------------------------------\n";
  Args args(argc, argv);
  args.ncrit = 4;
  args.numBodies = bodies.size();

  P = args.P;
  THETA = args.theta;
  PHI = args.phi;
  NCRIT = args.ncrit;
  VERBOSE = args.verbose;
  const int numBodies = args.numBodies;
  args.show();
  std::cout << "-------------------------------------\n";

  // Initializes the particles with 0 potential and force
  std::cout << "Initializing particles\n";
  initTarget(bodies);

  // Builds the tree
  std::cout << "Building tree\n";
  Cells cells = buildTree(bodies);

  // Calculates FMM
  std::cout << "Calculating FMM\n";
  initKernel();
  upwardPass(cells);
  horizontalPass(cells, cells);
  downwardPass(cells);

  // Searches for the new position of the probe
  std::cout << "Searching for probe index\n";
  int probe_i = -1;
  for(int i=0; i<np+1; ++i){
    if( bodies[i].X[0]==probe[0]
        && bodies[i].X[1]==probe[1]
        && bodies[i].X[2]==probe[2]){
      probe_i = i;
      break;
    }
  }
  if( probe_i==-1 ){
    std::cout << "\tCRITICAL: Probe not found!!!\n";
  }
  else{
    std::cout << "\tProbe found at i=" << probe_i << "\n";
  }


  std::cout << "\n**************************************\n";
  std::cout << "Potential at X=["<<bodies[probe_i].X[0]<<","<<bodies[probe_i].X[1]<<",";
  std::cout << bodies[probe_i].X[2] << "]\n";
  std::cout << "\tFMM potential\n\t";
  for(int ind=0; ind<3; ind++){
    std::cout << "\t" << bodies[probe_i].p[ind];
  }
  std::cout << "\n\tFMM potential converted into vortex potential\n\t";
  for(int ind=0; ind<3; ind++){
    std::cout << "\t" << (-1)/(4*M_PIl) * bodies[probe_i].p[ind];
  }
  std::cout << "\n**************************************\n";
  std::cout << "Jacobian at X=["<<bodies[probe_i].X[0]<<","<<bodies[probe_i].X[1]<<",";
  std::cout << bodies[probe_i].X[2] << "]\n";
  std::cout << "\tFMM Jacobian\n\t";
  for(int ind=0; ind<3; ind++){
    std::cout << "\t" << bodies[probe_i].J[ind*3 + 0];
    std::cout << "\t" << bodies[probe_i].J[ind*3 + 1];
    std::cout << "\t" << bodies[probe_i].J[ind*3 + 2] << "\n\t" ;
  }
  std::cout << "FMM Jacobian converted into vortex Jacobian\n\t";
  for(int ind=0; ind<3; ind++){
    std::cout << "\t" << (-1)/(4*M_PIl) * bodies[probe_i].J[ind*3 + 0];
    std::cout << "\t" << (-1)/(4*M_PIl) * bodies[probe_i].J[ind*3 + 1];
    std::cout << "\t" << (-1)/(4*M_PIl) * bodies[probe_i].J[ind*3 + 2] << "\n\t" ;
  }
  std::cout << "\n**************************************\n";

  return 0;
}

/* Calculates and returns the FMM's vector potential field and jacobian at
 each position of a particle field.

 INPUT:
  npX_n_param matrix `pfield`, where np is an arbitrary number of particles
  and _n_param is the number of parameters that defines a particle.

OUTPUT:
  npX(_n_param+3+9*4) matrix `OUT` where OUT[i][0:_n_param-1] is the input
  parameters of the i-th particle, OUT[i][_n_param:_n_param+2] is the vector
  potential at the position of the i-th particle, and
  OUT[i][_n_param+3:_n_param+3+9-1] is the linearized jacobian,
  OUT[i][_n_param+3+9:_n_param+3+9*2-1] is dJdx1,
  OUT[i][_n_param+3+9*2:_n_param+3+9*3-1] is dJdx2,
  OUT[i][_n_param+3+9*3:_n_param+3+9*4-1] is dJdx3.
  OUT[i][_n_param+3+9*4:_n_param+3+9*4+2] is the Particle Strength Exchange term.

NOTE: The particle order in `OUT` matches that of `pfield`.*/
void calculate(double** pfield, int np, double** OUT,
                double** OUT_lagrandist, int** OUT_lagrandist_np,
                int p=10, int ncrit=4, double theta=.4, double phi=.5, bool verbose=true,
                int p2p_type=0, int l2p_type=0, bool rbf=false, bool lagran_dist=false,
                double lambda=1.5, int max_lagrandist_np=0)
{
  // Dummy argument values
  char aux0 = 'a';
  char * aux1 = &aux0;
  int argc = 1;
  char ** argv = &aux1;

  // Generates bodies
  if(verbose) std::cout << "Generating bodies\n";
  Bodies bodies = genBodies(pfield, np);

  // Overwrites initialization values of ExaFMM
  if(verbose) std::cout << "Defining FMM parameters\n";
  if(verbose) std::cout << "-------------------------------------\n";
  // P = 10;
  // THETA = .4;
  // NCRIT = 4;
  P = p;
  THETA = theta;
  PHI = phi;
  NCRIT = ncrit;
  VERBOSE = verbose;
  P2P_TYPE = p2p_type;
  L2P_TYPE = l2p_type;
  RBF = rbf;
  LAMBDA = lambda;
  const int numBodies = bodies.size();
  if(verbose) std::cout << "\tP:\t"<<P<<"\n\tTheta:\t"<<THETA<<"\n\tPhi:\t"<<PHI;
  if(verbose) std::cout << "\n\tncrit:\t"<<NCRIT<<"\n\tnumBodies:\t"<<numBodies<<"\n";
  if(verbose) std::cout << "\n\tP2P_TYPE:\t"<<P2P_TYPE<<"\n\tL2P_TYPE:\t"<<L2P_TYPE<<"\n";
  if(verbose) std::cout << "-------------------------------------\n";

  // Initializes the particles with 0 potential and force
  if(verbose) std::cout << "Initializing particles\n";
  initTarget(bodies);

  // Builds the tree
  if(verbose) std::cout << "Building tree\n";
  Cells cells = buildTree(bodies);

  // Adds particles to take care of Lagrangian distortion
  int extra_np = 0;
  if(lagran_dist){

    // Generate extra particles
    Bodies lagran_bodies = closeConnectivity(cells);

    for(size_t i=0; i<lagran_bodies.size(); i++){
      // Gives them an index
      lagran_bodies[i].index = np+i+1;
      // Makes sure they have no strength
      lagran_bodies[i].q = 0;
    }

    // Adds them to the initial ones
    bodies.insert(bodies.end(), lagran_bodies.begin(), lagran_bodies.end());

    // Regenerates the tree
    cells = buildTree(bodies);

    extra_np = lagran_bodies.size();
  }

  if(extra_np>max_lagrandist_np) std::cout <<
                    "ERROR: Maximum number of Lagrangian distortion particles "
                    << max_lagrandist_np << " reached!" << "\n";

  // Calculates FMM
  if(verbose) std::cout << "Calculating FMM\n";
  initKernel();
  upwardPass(cells);
  // clock_t begin = std::clock();
  horizontalPass(cells, cells);
  // clock_t end = std::clock();
  // std::cout << "Time:"<<double(end-begin)/ CLOCKS_PER_SEC<<"\n";
  downwardPass(cells);


  // Prepares output
  for(int i=0; i<np+extra_np; i++){
    int index = bodies[i].index-1;
    if(index < np){
      OUT[index][_x_i-1] = bodies[i].X[0];
      OUT[index][_y_i-1] = bodies[i].X[1];
      OUT[index][_z_i-1] = bodies[i].X[2];
      OUT[index][_gamma_x_i-1] = bodies[i].q[0];
      OUT[index][_gamma_y_i-1] = bodies[i].q[1];
      OUT[index][_gamma_z_i-1] = bodies[i].q[2];
      OUT[index][_blob_sigma_i-1] = bodies[i].sigma;
      OUT[index][_vol_i-1] = bodies[i].vol;
      for(int ind=0; ind<3; ind++) OUT[index][ (_n_param) + ind] = bodies[i].p[ind];
      for(int ind=0; ind<9; ind++) OUT[index][ (_n_param+3) + ind] = bodies[i].J[ind];
      for(int ind=0; ind<9; ind++) OUT[index][ (_n_param+3+9) + ind] = bodies[i].dJdx1[ind];
      for(int ind=0; ind<9; ind++) OUT[index][ (_n_param+3+9*2) + ind] = bodies[i].dJdx2[ind];
      for(int ind=0; ind<9; ind++) OUT[index][ (_n_param+3+9*3) + ind] = bodies[i].dJdx3[ind];
      for(int ind=0; ind<3; ind++) OUT[index][ (_n_param+3+9*4) + ind] = bodies[i].pse[ind];
    }
    else{
      index = index-np;
      OUT_lagrandist[index][_x_i-1] = bodies[i].X[0];
      OUT_lagrandist[index][_y_i-1] = bodies[i].X[1];
      OUT_lagrandist[index][_z_i-1] = bodies[i].X[2];
      OUT_lagrandist[index][_gamma_x_i-1] = bodies[i].q[0];
      OUT_lagrandist[index][_gamma_y_i-1] = bodies[i].q[1];
      OUT_lagrandist[index][_gamma_z_i-1] = bodies[i].q[2];
      OUT_lagrandist[index][_blob_sigma_i-1] = bodies[i].sigma;
      OUT_lagrandist[index][_vol_i-1] = bodies[i].vol;
      for(int ind=0; ind<3; ind++) OUT_lagrandist[index][ (_n_param) + ind] = bodies[i].p[ind];
      for(int ind=0; ind<9; ind++) OUT_lagrandist[index][ (_n_param+3) + ind] = bodies[i].J[ind];
      for(int ind=0; ind<9; ind++) OUT_lagrandist[index][ (_n_param+3+9) + ind] = bodies[i].dJdx1[ind];
      for(int ind=0; ind<9; ind++) OUT_lagrandist[index][ (_n_param+3+9*2) + ind] = bodies[i].dJdx2[ind];
      for(int ind=0; ind<9; ind++) OUT_lagrandist[index][ (_n_param+3+9*3) + ind] = bodies[i].dJdx3[ind];
      for(int ind=0; ind<3; ind++) OUT_lagrandist[index][ (_n_param+3+9*4) + ind] = bodies[i].pse[ind];
    }
  }
  if(lagran_dist) OUT_lagrandist_np[0][0] = extra_np;
}


void testcalculate(double** pfield, int np, double** OUT)
{
  for(int i=0; i<np; i++){
    OUT[i][0] = i;
    for(int j=0; j<_n_param; j++){
      OUT[i][j+1] = pfield[i][j];
    }
  }
}

void testArray(int** array, unsigned height, unsigned width)
{
  int i = 0;
  for (int h = 0; h < height; h++)
  {
    for (int w = 0; w < width; w++)
    {
      array[h][w] = i++;
    }
  }
}


// CxxWrap test function
std::string greet()
{
   return "hello, world";
}

// Makes available the module CppExaFMM in julia
JULIA_CPP_MODULE_BEGIN(registry)
  jlcxx::Module& fmm = registry.create_module("CppExaFMM");
  fmm.method("greet", &greet);
  fmm.method("main", &main);
  fmm.method("testArray", &testArray);
  fmm.method("testcalculate", &testcalculate);
  fmm.method("calculate", &calculate);
JULIA_CPP_MODULE_END
