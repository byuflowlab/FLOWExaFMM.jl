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

    mod.method("greet", &greet);
}
