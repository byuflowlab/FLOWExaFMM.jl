/**
 * @Author: Eduardo Alvarez <user>
 * @Date:   2020-08-1T16:10:58-06:00
 * @Email:  Edo.AlvarezR@gmail.com
 * @Comments: I replaced fmm.cxx from the original build to make this module
 *            an interface between FLOWVPM and exafmm. The original module can be
 *            found at fmm-org.cxx
 */

 #include "jlcxx/jlcxx.hpp"  //C++ wrapper for julia

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

// Extra tools for interacting
#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>  //for 'typeid' to work
#include <string>

#define M_PIl          3.141592653589793238462643383279502884L


//! Initialize bodies with zero/null values
void initBodiesZero(Bodies & bodies) {
  for (size_t b=0; b!=bodies.size(); ++b) {
    bodies[b].X = 0;
    bodies[b].q = 0;
    bodies[b].vort = 0;
    bodies[b].sigma = -1;
    bodies[b].vol = -1;
    bodies[b].index = -1;
    bodies[b].p = 0;
    bodies[b].J = 0;
    bodies[b].dJdx1 = 0;
    bodies[b].dJdx2 = 0;
    bodies[b].dJdx3 = 0;
    bodies[b].pse = 0;
  }
}

//! Create and initialize an array of nb bodies
Bodies genBodies(int nb){

  Bodies bodies(nb);
  initBodiesZero(bodies);

  return bodies;
}

//! Return the i-th body (0-indexed)
Body* getBody(Bodies & bodies, int i){return &bodies[i];}

// Test function
std::string greet()
{
   return "hello, world";
}

// ! Overwrite the target body with the source body
void overwriteBody(Bodies & bodies, int trg, int src){
  bodies[trg].X = bodies[src].X;
  bodies[trg].q = bodies[src].q;
  bodies[trg].sigma = bodies[src].sigma;
  bodies[trg].vol = bodies[src].vol;
  bodies[trg].index = bodies[src].index;
}

bool getPrecision(){
  return single_prec;
}

void calculate(Bodies & bodies, int np, int p, int ncrit, real_t theta, real_t phi,
                bool verbose, int p2p_type, int l2p_type, bool rbf){

    // Dummy argument values
    char aux0 = 'a';
    char * aux1 = &aux0;
    int argc = 1;
    char ** argv = &aux1;

    // Overwrite initialization values of ExaFMM
    if(verbose) std::cout << "Defining FMM parameters\n";
    if(verbose) std::cout << "-------------------------------------\n";
    P = p;
    THETA = theta;
    PHI = phi;
    NCRIT = ncrit;
    VERBOSE = verbose;
    P2P_TYPE = p2p_type;
    L2P_TYPE = l2p_type;
    RBF = rbf;
    const int numBodies = np;
    if(verbose) std::cout << "\tP:\t"<<P<<"\n\tTheta:\t"<<THETA<<"\n\tPhi:\t"<<PHI;
    if(verbose) std::cout << "\n\tncrit:\t"<<NCRIT<<"\n\tnumBodies:\t"<<numBodies<<"\n";
    if(verbose) std::cout << "\n\tP2P_TYPE:\t"<<P2P_TYPE<<"\n\tL2P_TYPE:\t"<<L2P_TYPE<<"\n";
    if(verbose) std::cout << "-------------------------------------\n";

    if(verbose){
      std::cout << "-------------------------------------\n";
      for(int i=0; i<np; i++){
          std::cout << "Particle #" << i+1 <<"\n";
          std::cout << "\tindex = " << bodies[i].index[0] << "\n";
          std::cout << "\tsigma = " << bodies[i].sigma[0] << "\n";
          std::cout << "\tX = [" << bodies[i].X[0]<<", "<<bodies[i].X[1]<<", "<<bodies[i].X[2]<<"]\n";
          std::cout << "\tq = [" << bodies[i].q[0]<<", "<<bodies[i].q[1]<<", "<<bodies[i].q[2]<<"]\n";
          std::cout << "\tJ = [";
          for(int j=0; j<9; j++){
              std::cout << bodies[i].J[j];
              if(j<8) std::cout << ", ";
              else std::cout << "]\n";
          }
      }
      std::cout << "-------------------------------------\n";
    }

    // Initialize the particles with 0 potential and force
    if(verbose) std::cout << "Initializing particles\n";
    initTarget(bodies, np);

    // Build the tree
    if(verbose) std::cout << "Building tree\n";
    Cells cells = buildTree(bodies, np);

    // Calculate FMM
    if(verbose) std::cout << "Calculating FMM\n";

    initKernel();
    upwardPass(cells);
    // clock_t begin = std::clock();
    horizontalPass(cells, cells);
    // clock_t end = std::clock();
    // std::cout << "Time:"<<double(end-begin)/ CLOCKS_PER_SEC<<"\n";
    downwardPass(cells);

    if(verbose){
      std::cout << "-------------------------------------\n";
      for(int i=0; i<np; i++){
          std::cout << "Particle #" << i+1 <<"\n";
          std::cout << "\tindex = " << bodies[i].index[0] << "\n";
          std::cout << "\tX = [" << bodies[i].X[0]<<", "<<bodies[i].X[1]<<", "<<bodies[i].X[2]<<"]\n";
          std::cout << "\tq = [" << bodies[i].q[0]<<", "<<bodies[i].q[1]<<", "<<bodies[i].q[2]<<"]\n";
          std::cout << "\tJ = [";
          for(int j=0; j<9; j++){
              std::cout << bodies[i].J[j];
              if(j<8) std::cout << ", ";
              else std::cout << "]\n";
          }
      }
      std::cout << "-------------------------------------\n";
    }

    // Calculate FMM
    if(verbose) std::cout << "Done!\n";
}

// Exposing types and functions to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<Body>("Body")
      .constructor()

      .method("get_Xref", &Body::get_Xref)
      .method("get_Xi", &Body::get_Xi)
      .method("set_Xi", &Body::set_Xi)

      .method("get_qref", &Body::get_qref)
      .method("get_qi", &Body::get_qi)
      .method("set_qi", &Body::set_qi)

      .method("get_vortref", &Body::get_vortref)
      .method("get_vorti", &Body::get_vorti)
      .method("set_vorti", &Body::set_vorti)

      .method("get_sigmaref", &Body::get_sigmaref)
      .method("get_sigmai", &Body::get_sigmai)
      .method("set_sigmai", &Body::set_sigmai)

      .method("get_volref", &Body::get_volref)
      .method("get_voli", &Body::get_voli)
      .method("set_voli", &Body::set_voli)

      .method("get_indexref", &Body::get_indexref)
      .method("get_indexi", &Body::get_indexi)
      .method("set_indexi", &Body::set_indexi)

      .method("get_pref", &Body::get_pref)
      .method("get_pi", &Body::get_pi)
      .method("set_pi", &Body::set_pi)

      .method("get_Jref", &Body::get_Jref)
      .method("get_Ji", &Body::get_Ji)
      .method("set_Ji", &Body::set_Ji)

      .method("get_dJdx1ref", &Body::get_dJdx1ref)
      .method("get_dJdx1i", &Body::get_dJdx1i)
      .method("set_dJdx1i", &Body::set_dJdx1i)

      .method("get_dJdx2ref", &Body::get_dJdx2ref)
      .method("get_dJdx2i", &Body::get_dJdx2i)
      .method("set_dJdx2i", &Body::set_dJdx2i)

      .method("get_dJdx3ref", &Body::get_dJdx3ref)
      .method("get_dJdx3i", &Body::get_dJdx3i)
      .method("set_dJdx3i", &Body::set_dJdx3i)

      .method("get_pseref", &Body::get_pseref)
      .method("get_psei", &Body::get_psei)
      .method("set_psei", &Body::set_psei);

    mod.add_type<Bodies>("Bodies");
    mod.method("genBodies", &genBodies);
    mod.method("getBody", &getBody);
    mod.method("overwriteBody", &overwriteBody);

    mod.method("getPrecision", &getPrecision);

    mod.method("calculate", &calculate);

    mod.method("greet", &greet);
}
