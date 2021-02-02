/**
 * @Author: Eduardo Alvarez <user>
 * @Date:   2020-08-1T16:10:58-06:00
 * @Email:  Edo.AlvarezR@gmail.com
 * @Comments: I replaced fmm.cxx from the original build to make this module
 *            an interface between FLOWVPM and exafmm. The original module can be
 *            found at fmm-org.cxx
 */

 #include "jlcxx/jlcxx.hpp"  //C++ wrapper for julia


 // Extra tools for interacting
 #include <iostream>
 #include <fstream>
 #include <sstream>
 #include <string>
 #include <typeinfo>  //for 'typeid' to work


 // ExaFMM modules
#include "exafmm.h"
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
  // bodies[trg].vort = bodies[src].vort;
  bodies[trg].sigma = bodies[src].sigma;
  bodies[trg].vol = bodies[src].vol;
  bodies[trg].index = bodies[src].index;
  // bodies[trg].p = bodies[src].p;
  // bodies[trg].J = bodies[src].J;
  // bodies[trg].dJdx1 = bodies[src].dJdx1;
  // bodies[trg].dJdx2 = bodies[src].dJdx2;
  // bodies[trg].dJdx3 = bodies[src].dJdx3;
  // bodies[trg].pse = bodies[src].pse;
}

void sortBodies(Bodies & bodies, int nb){
  vec3 X;                                     //!< Position
  vec3 q;                                     //!< Vector charge
  vec3 vort;                                  //!< Vorticity for RBF
  vec1 sigma;                                 //!< Vortex blob radius
  vec1 vol;                                   //!< Vortex particle volume (it will also use this for RBF)
  vec1int index;                              //!< Index of this body
  vec3 p;                                     //!< Vector potential
  vec9 J;                                     //!< Jacobian
  vec9 dJdx1;                                 //!< Partial of Jacobian to x1
  vec9 dJdx2;                                 //!< Partial of Jacobian to x2
  vec9 dJdx3;                                 //!< Partial of Jacobian to x3
  vec3 pse;                                   //!< Particle Strength Exchange
  int ind;

  for(int i=0; i<nb; i++){ // Iterate over array elements
    while(bodies[i].index[0] != i){ // If position has the wrong element

      // Lift up element in this position
      X = bodies[i].X;
      q = bodies[i].q;
      vort = bodies[i].vort;
      sigma = bodies[i].sigma;
      vol = bodies[i].vol;
      index = bodies[i].index;
      p = bodies[i].p;
      J = bodies[i].J;
      dJdx1 = bodies[i].dJdx1;
      dJdx2 = bodies[i].dJdx2;
      dJdx3 = bodies[i].dJdx3;
      pse = bodies[i].pse;

      ind = bodies[i].index[0];

      // Move element that is in the position where this element should go to this position
      bodies[i].X = bodies[ind].X;
      bodies[i].q = bodies[ind].q;
      bodies[i].vort = bodies[ind].vort;
      bodies[i].sigma = bodies[ind].sigma;
      bodies[i].vol = bodies[ind].vol;
      bodies[i].index = bodies[ind].index;
      bodies[i].p = bodies[ind].p;
      bodies[i].J = bodies[ind].J;
      bodies[i].dJdx1 = bodies[ind].dJdx1;
      bodies[i].dJdx2 = bodies[ind].dJdx2;
      bodies[i].dJdx3 = bodies[ind].dJdx3;
      bodies[i].pse = bodies[ind].pse;

      // Place this element where it should go
      bodies[ind].X = X;
      bodies[ind].q = q;
      bodies[ind].vort = vort;
      bodies[ind].sigma = sigma;
      bodies[ind].vol = vol;
      bodies[ind].index = index;
      bodies[ind].p = p;
      bodies[ind].J = J;
      bodies[ind].dJdx1 = dJdx1;
      bodies[ind].dJdx2 = dJdx2;
      bodies[ind].dJdx3 = dJdx3;
      bodies[ind].pse = pse;
    }
  }
}

bool getPrecision(){
  return single_prec;
}

void calculate(Bodies & bodies, int np, int p, int ncrit, real_t theta, real_t phi,
                bool verbose, int p2p_type, int l2p_type,
                bool rbf, bool sgs, bool transposed, bool reset, bool reset_sgs,
                bool sort=true){

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
    SGS = sgs;
    TRANSPOSED = transposed;
    const int numBodies = np;
    if(verbose) std::cout << "\tP:\t"<<P<<"\n\tTheta:\t"<<THETA<<"\n\tPhi:\t"<<PHI;
    if(verbose) std::cout << "\n\tncrit:\t"<<NCRIT<<"\n\tnumBodies:\t"<<numBodies<<"\n";
    if(verbose) std::cout << "\n\tP2P_TYPE:\t"<<P2P_TYPE<<"\n\tL2P_TYPE:\t"<<L2P_TYPE<<"\n";
    if(verbose) std::cout << "-------------------------------------\n";

    // Index particles according to array position
    if(sort){
      for(int i=0; i<np; i++){
        bodies[i].index[0] = i;
      }
    }

    if(verbose){
      std::cout << "-------------------------------------\n";
      // std::cout << std::fixed;
      std::cout << std::setprecision(17);
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
    if(reset){
      if(verbose) std::cout << "Initializing particles\n";
      initTarget(bodies, np);
    }
    if(reset_sgs){
      if(verbose) std::cout << "Initializing particles SGS\n";
      initTarget2(bodies, np);
    }

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

    // NOTE: If the bodies are not sorted back to their original position in the
    //  memory, somehow the single vortex ring case falls apart due to small
    //  floating-point-precision numerical randomness.
    if(sort) sortBodies(bodies, np);

    if(verbose){
      std::cout << "-------------------------------------\n";
      // std::cout << std::fixed;
      std::cout << std::setprecision(17);
      for(int i=0; i<np; i++){
          std::cout << "Particle #" << i+1 <<"\n";
          std::cout << "\tindex = " << bodies[i].index[0] << "\n";
          std::cout << "\tX = [" << bodies[i].X[0]<<", "<<bodies[i].X[1]<<", "<<bodies[i].X[2]<<"]\n";
          std::cout << "\tq = [" << bodies[i].q[0]<<", "<<bodies[i].q[1]<<", "<<bodies[i].q[2]<<"]\n";
          std::cout << "\ts = " << bodies[i].sigma[0] << "\n";
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
