#ifndef rbf_h
#define rbf_h
#include "exafmm.h"

namespace exafmm {



  /* Radial basis function approximation

    INPUTS
    * bodies:     Particles.
    * cells:      Tree structure with all the particles.
    * vortfield:  npx3 matrix with the vorticity at each particle, with
                  vortfield[body.index-1] being the vorticity at the position
                  of particle of index `body.index`.
    * np:         Number of particles.
    * itmax:      Maximum number of iterations.
    * tol:        Tolerance.
  */
  void rbf(Bodies & bodies, Cells & cells, real_t** vortfield, int np,
                                                        int itmax, real_t tol) {
    // const int itmax = 5;
    // const float tol = 1e-3;

    // Initial guess: vorticity*vol
    for(int i=0; i<np; i++){
      for(int j=0; j<3; j++) bodies[i].q[j] = vortfield[i][j]*bodies[i].vol
    }


  }

}
#endif
