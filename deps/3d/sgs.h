
#ifndef sgs_h
#define sgs_h
#include "exafmm.h"

namespace exafmm {

    // Vortex strethcing subgrid-scale model on Winckelman's algebraic kernel (See notebook 20210201)
    void SGS_stretch1_P2P_alg(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          // Backscatter filter criterion f_pq
          if ( (Bi[i].q[0]*S[0] + Bi[i].q[1]*S[1] + Bi[i].q[2]*S[2]) >= 0) {

            real_t R2 = norm(Bi[i].X - Bj[j].X);
            real_t aux2 = R2 / (Bj[j].sigma[0]*Bj[j].sigma[0]);
            real_t zeta = 7.5 / std::pow(aux2 + 1, 3.5) / (Bj[j].sigma[0]*Bj[j].sigma[0]*Bj[j].sigma[0]) / (M_PI*4);

            // Accumulate zeta * S
            // (Divide S by 4*pi to obtain velocity derivatives)
            SGS[0] += aux1*zeta*S[0];
            SGS[1] += aux1*zeta*S[1];
            SGS[2] += aux1*zeta*S[2];
          }

        }

        // Save it where FLOWVPM will find it
        Bi[i].J[0] -= SGS[0];
        Bi[i].J[1] -= SGS[1];
        Bi[i].J[2] -= SGS[2];
      }
    }

    // Vortex strethcing subgrid-scale model on Gaussian-Erf kernel (See notebook 20210201)
    void SGS_stretch1_P2P_Gaussian(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          // Backscatter filter criterion f_pq
          if (Bi[i].q[0]*S[0] + Bi[i].q[1]*S[1] + Bi[i].q[2]*S[2] >= 0) {

            real_t R2 = norm(Bi[i].X - Bj[j].X);
            real_t S2 = 2 * Bj[j].sigma[0] * Bj[j].sigma[0];
            real_t zeta = std::exp(-R2 / S2) / (M_PI * S2) / sqrt(M_PI * S2);

            // Accumulate zeta * S
            // (Divide S by 4*pi to obtain velocity derivatives)
            SGS[0] += aux1*zeta*S[0];
            SGS[1] += aux1*zeta*S[1];
            SGS[2] += aux1*zeta*S[2];
          }

        }

        // Save it where FLOWVPM will find it
        Bi[i].J[0] -= SGS[0];
        Bi[i].J[1] -= SGS[1];
        Bi[i].J[2] -= SGS[2];
      }
    }


    // Vortex strethcing subgrid-scale model on Winckelman's algebraic kernel (See notebook 20210201)
    void SGS_stretch2_P2P_alg(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          real_t R2 = norm(Bi[i].X - Bj[j].X);
          real_t aux2 = R2 / (Bj[j].sigma[0]*Bj[j].sigma[0]);
          real_t zeta = 7.5 / std::pow(aux2 + 1, 3.5) / (Bj[j].sigma[0]*Bj[j].sigma[0]*Bj[j].sigma[0]) / (M_PI*4);

          // Accumulate zeta * S
          // (Divide S by 4*pi to obtain velocity derivatives)
          SGS[0] += aux1*zeta*S[0];
          SGS[1] += aux1*zeta*S[1];
          SGS[2] += aux1*zeta*S[2];

        }

        // Backscatter filter criterion f_pq
        if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] -= SGS[0];
          Bi[i].J[1] -= SGS[1];
          Bi[i].J[2] -= SGS[2];
        }
      }
    }

    // Vortex strethcing subgrid-scale model on Gaussian-Erf kernel (See notebook 20210201)
    void SGS_stretch2_P2P_Gaussian(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          real_t R2 = norm(Bi[i].X - Bj[j].X);
          real_t S2 = 2 * Bj[j].sigma[0] * Bj[j].sigma[0];
          real_t zeta = std::exp(-R2 / S2) / (M_PI * S2) / sqrt(M_PI * S2);

          // Accumulate zeta * S
          // (Divide S by 4*pi to obtain velocity derivatives)
          SGS[0] += aux1*zeta*S[0];
          SGS[1] += aux1*zeta*S[1];
          SGS[2] += aux1*zeta*S[2];

        }

        // Backscatter filter criterion f_pq
        if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] -= SGS[0];
          Bi[i].J[1] -= SGS[1];
          Bi[i].J[2] -= SGS[2];
        }
      }
    }


    // Vortex strethcing subgrid-scale model on Winckelman's algebraic kernel (See notebook 20210201)
    void SGS_stretch0_P2P_alg(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          real_t R2 = norm(Bi[i].X - Bj[j].X);
          real_t aux2 = R2 / (Bj[j].sigma[0]*Bj[j].sigma[0]);
          real_t zeta = 7.5 / std::pow(aux2 + 1, 3.5) / (Bj[j].sigma[0]*Bj[j].sigma[0]*Bj[j].sigma[0]) / (M_PI*4);

          // Accumulate zeta * S
          // (Divide S by 4*pi to obtain velocity derivatives)
          SGS[0] += aux1*zeta*S[0];
          SGS[1] += aux1*zeta*S[1];
          SGS[2] += aux1*zeta*S[2];

        }

        // Backscatter filter criterion f_pq
        // if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] -= SGS[0];
          Bi[i].J[1] -= SGS[1];
          Bi[i].J[2] -= SGS[2];
        // }
      }
    }

    // Vortex strethcing subgrid-scale model on Gaussian-Erf kernel (See notebook 20210201)
    void SGS_stretch0_P2P_Gaussian(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
          // NOTE: this is negative since U is flipped in the VPM
          if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }
          else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
            // x-component
            S[0] = -Bj[j].q[0]*( (Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2]) - (Bj[j].dJdx1[2*3 + 1]-Bj[j].dJdx1[1*3 + 2]) );
            S[0] -= Bj[j].q[1]*( (Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2]) - (Bj[j].dJdx2[2*3 + 1]-Bj[j].dJdx2[1*3 + 2]) );
            S[0] -= Bj[j].q[2]*( (Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2]) - (Bj[j].dJdx3[2*3 + 1]-Bj[j].dJdx3[1*3 + 2]) );
            // y-component
            S[1] = -Bj[j].q[0]*( (Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0]) - (Bj[j].dJdx1[0*3 + 2]-Bj[j].dJdx1[2*3 + 0]) );
            S[1] -= Bj[j].q[1]*( (Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0]) - (Bj[j].dJdx2[0*3 + 2]-Bj[j].dJdx2[2*3 + 0]) );
            S[1] -= Bj[j].q[2]*( (Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0]) - (Bj[j].dJdx3[0*3 + 2]-Bj[j].dJdx3[2*3 + 0]) );
            // z-component
            S[2] = -Bj[j].q[0]*( (Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1]) - (Bj[j].dJdx1[1*3 + 0]-Bj[j].dJdx1[0*3 + 1]) );
            S[2] -= Bj[j].q[1]*( (Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1]) - (Bj[j].dJdx2[1*3 + 0]-Bj[j].dJdx2[0*3 + 1]) );
            S[2] -= Bj[j].q[2]*( (Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1]) - (Bj[j].dJdx3[1*3 + 0]-Bj[j].dJdx3[0*3 + 1]) );
          }

          // I am not sure why the calculation doesn't match the FLOWVPM's
          // calculation unless I multiply it by -1 here
          S *= -1;

          real_t R2 = norm(Bi[i].X - Bj[j].X);
          real_t S2 = 2 * Bj[j].sigma[0] * Bj[j].sigma[0];
          real_t zeta = std::exp(-R2 / S2) / (M_PI * S2) / sqrt(M_PI * S2);

          // Accumulate zeta * S
          // (Divide S by 4*pi to obtain velocity derivatives)
          SGS[0] += aux1*zeta*S[0];
          SGS[1] += aux1*zeta*S[1];
          SGS[2] += aux1*zeta*S[2];

        }

        // Backscatter filter criterion f_pq
        // if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] -= SGS[0];
          Bi[i].J[1] -= SGS[1];
          Bi[i].J[2] -= SGS[2];
        // }
      }
    }


    // Vortex strethcing subgrid-scale model on Winckelman's algebraic kernel (See notebook 20210201)
    void SGS_M2_P2P_alg(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          real_t R2 = norm(Bi[i].X - Bj[j].X);

          if(R2 != 0){

            // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
            // NOTE: this is negative since U is flipped in the VPM
            if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
              // x-component
              S[0] = -Bj[j].q[0]*( Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2] );
              S[0] -= Bj[j].q[1]*( Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0] );
              S[0] -= Bj[j].q[2]*( Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1] );
              // y-component
              S[1] = -Bj[j].q[0]*( Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2] );
              S[1] -= Bj[j].q[1]*( Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0] );
              S[1] -= Bj[j].q[2]*( Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1] );
              // z-component
              S[2] = -Bj[j].q[0]*( Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2] );
              S[2] -= Bj[j].q[1]*( Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0] );
              S[2] -= Bj[j].q[2]*( Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1] );
            }
            else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
              // x-component
              S[0] = -Bj[j].q[0]*( Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2] );
              S[0] -= Bj[j].q[1]*( Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2] );
              S[0] -= Bj[j].q[2]*( Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2] );
              // y-component
              S[1] = -Bj[j].q[0]*( Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0] );
              S[1] -= Bj[j].q[1]*( Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0] );
              S[1] -= Bj[j].q[2]*( Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0] );
              // z-component
              S[2] = -Bj[j].q[0]*( Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1] );
              S[2] -= Bj[j].q[1]*( Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1] );
              S[2] -= Bj[j].q[2]*( Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1] );
            }

            // I am not sure why the calculation doesn't match the FLOWVPM's
            // calculation unless I multiply it by -1 here
            S *= -1;

            real_t aux2 = R2 / (Bj[j].sigma[0]*Bj[j].sigma[0]);
            real_t zeta = 7.5 / std::pow(aux2 + 1, 3.5) / (Bj[j].sigma[0]*Bj[j].sigma[0]*Bj[j].sigma[0]) / (M_PI*4);

            // Accumulate zeta * S
            // (Divide S by 4*pi to obtain velocity derivatives)
            SGS[0] += aux1*zeta*S[0];
            SGS[1] += aux1*zeta*S[1];
            SGS[2] += aux1*zeta*S[2];

          }

        }

        // Backscatter filter criterion f_pq
        // if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] += SGS[0];
          Bi[i].J[1] += SGS[1];
          Bi[i].J[2] += SGS[2];
        // }
      }
    }

    // Vortex strethcing subgrid-scale model on Gaussian-Erf kernel (See notebook 20210201)
    void SGS_M2_P2P_Gaussian(Cell * Ci, Cell * Cj) {
      Body * Bi = Ci->body;
      Body * Bj = Cj->body;

      vec3 S = 0;   // Store stretching here
      real_t aux1 = 1/(4*M_PI);

      // Iterate over p particle
      for (int i=0; i<Ci->numBodies; i++) {
        vec3 SGS = 0;

        // Influence of q particle on p  (or j on i)
        for (int j=0; j<Cj->numBodies; j++) {

          real_t R2 = norm(Bi[i].X - Bj[j].X);

          if(R2 != 0){

            // Calculate stretching between p and q as S = (Gamma_q dot nabla) (u_p - x_q) )
            // NOTE: this is negative since U is flipped in the VPM
            if (TRANSPOSED){         // S = (Gamma_q dot nabla') (u_p - x_q) )
              // x-component
              S[0] = -Bj[j].q[0]*( Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2] );
              S[0] -= Bj[j].q[1]*( Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0] );
              S[0] -= Bj[j].q[2]*( Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1] );
              // y-component
              S[1] = -Bj[j].q[0]*( Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2] );
              S[1] -= Bj[j].q[1]*( Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0] );
              S[1] -= Bj[j].q[2]*( Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1] );
              // z-component
              S[2] = -Bj[j].q[0]*( Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2] );
              S[2] -= Bj[j].q[1]*( Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0] );
              S[2] -= Bj[j].q[2]*( Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1] );
            }
            else{                   // S = (Gamma_q dot nabla) (u_p - x_q) )
              // x-component
              S[0] = -Bj[j].q[0]*( Bi[i].dJdx1[2*3 + 1]-Bi[i].dJdx1[1*3 + 2] );
              S[0] -= Bj[j].q[1]*( Bi[i].dJdx2[2*3 + 1]-Bi[i].dJdx2[1*3 + 2] );
              S[0] -= Bj[j].q[2]*( Bi[i].dJdx3[2*3 + 1]-Bi[i].dJdx3[1*3 + 2] );
              // y-component
              S[1] = -Bj[j].q[0]*( Bi[i].dJdx1[0*3 + 2]-Bi[i].dJdx1[2*3 + 0] );
              S[1] -= Bj[j].q[1]*( Bi[i].dJdx2[0*3 + 2]-Bi[i].dJdx2[2*3 + 0] );
              S[1] -= Bj[j].q[2]*( Bi[i].dJdx3[0*3 + 2]-Bi[i].dJdx3[2*3 + 0] );
              // z-component
              S[2] = -Bj[j].q[0]*( Bi[i].dJdx1[1*3 + 0]-Bi[i].dJdx1[0*3 + 1] );
              S[2] -= Bj[j].q[1]*( Bi[i].dJdx2[1*3 + 0]-Bi[i].dJdx2[0*3 + 1] );
              S[2] -= Bj[j].q[2]*( Bi[i].dJdx3[1*3 + 0]-Bi[i].dJdx3[0*3 + 1] );
            }

            // I am not sure why the calculation doesn't match the FLOWVPM's
            // calculation unless I multiply it by -1 here
            S *= -1;

            real_t S2 = 2 * Bj[j].sigma[0] * Bj[j].sigma[0];
            real_t zeta = std::exp(-R2 / S2) / (M_PI * S2) / sqrt(M_PI * S2);

            // Accumulate zeta * S
            // (Divide S by 4*pi to obtain velocity derivatives)
            SGS[0] += aux1*zeta*S[0];
            SGS[1] += aux1*zeta*S[1];
            SGS[2] += aux1*zeta*S[2];
          }

        }

        // Backscatter filter criterion f_pq
        // if ( (Bi[i].q[0]*SGS[0] + Bi[i].q[1]*SGS[1] + Bi[i].q[2]*SGS[2]) >= 0) {
          // Save it where FLOWVPM will find it
          Bi[i].J[0] += SGS[0];
          Bi[i].J[1] += SGS[1];
          Bi[i].J[2] += SGS[2];
        // }
      }
    }

}
#endif
