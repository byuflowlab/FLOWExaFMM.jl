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

// Exposing types and functions to Julia
JLCXX_MODULE define_julia_module(jlcxx::Module& mod)
{
    mod.add_type<Body>("Body")
      .constructor()
      .method("get_Xref", &Body::get_Xref)
      .method("get_qref", &Body::get_qref)
      .method("get_vortref", &Body::get_vortref)
      .method("get_sigma", &Body::get_sigma)
      .method("get_vol", &Body::get_vol)
      .method("get_index", &Body::get_index)
      .method("get_pref", &Body::get_pref)
      .method("get_Jref", &Body::get_Jref)
      .method("get_dJdx1ref", &Body::get_dJdx1ref)
      .method("get_dJdx2ref", &Body::get_dJdx2ref)
      .method("get_dJdx3ref", &Body::get_dJdx3ref)
      .method("get_pseref", &Body::get_pseref)
      .method("set_sigma", &Body::set_sigma)
      .method("set_vol", &Body::set_vol)
      .method("set_index", &Body::set_index);

    mod.add_type<Bodies>("Bodies");
    mod.method("initBodiesZero", &initBodiesZero);
    mod.method("genBodies", &genBodies);
    mod.method("getBody", &getBody);
    mod.method("greet", &greet);
}
