// REFERENCES
// 1. Yokota, 2010, *Comparing the treecode with FMM on GPUs for vortex particles
//    simulations of a leapfrogging vortex ring*.

#ifndef kernel_h
#define kernel_h
#include "exafmm.h"

namespace exafmm {
  const complex_t I(0.,1.);

  inline int oddOrEven(int n) {
    return (((n) & 1) == 1) ? -1 : 1;
  }

  inline int ipow2n(int n) {
    return (n >= 0) ? 1 : oddOrEven(n);
  }

  //! Get r,theta,phi from x,y,z
  void cart2sph(const vec3 & dX, real_t & r, real_t & theta, real_t & phi) {
    r = sqrt(norm(dX));                                         // r = sqrt(x^2 + y^2 + z^2)
    theta = r<=EPS ? 0 : acos(dX[2] / r);                       // theta = acos(z / r)
    phi = r<=EPS ? 0 : atan2(dX[1], dX[0]);                     // phi = atan(y / x)
  }
  void cart2sph(const cvec3 & dX, cplx & r, cplx & theta, cplx & phi) {
    r = sqrt(norm(dX));
    theta = real_t(r)<=EPS ? cplx(0,0) : acos(dX[2] / r);
    phi = real_t(r)<=EPS ? cplx(0,0) : atan2(dX[1], dX[0]);
  }

  //! Spherical to cartesian coordinates
  void sph2cart(real_t r, real_t theta, real_t phi, const vec3 & spherical, vec3 & cartesian) {
    if (r<=EPS){
      std::cout << "Logic Error: encountered spherical to cartesian transformation with r==0!"<<
                      " This is most likely due to L2P with P in center of L"<<"\n";
      throw std::invalid_argument( "Invalid r==0" );
    }
    cartesian[0] = std::sin(theta) * std::cos(phi) * spherical[0] // x component (not x itself)
      + std::cos(theta) * std::cos(phi) / r * spherical[1]
      - std::sin(phi) / r / std::sin(theta) * spherical[2];
    cartesian[1] = std::sin(theta) * std::sin(phi) * spherical[0] // y component (not y itself)
      + std::cos(theta) * std::sin(phi) / r * spherical[1]
      + std::cos(phi) / r / std::sin(theta) * spherical[2];
    cartesian[2] = std::cos(theta) * spherical[0]                 // z component (not z itself)
      - std::sin(theta) / r * spherical[1];
  }

  //! Spherical to cartesian coordinates for CSDA
  void sph2cart(cplx r, cplx theta, cplx phi, const cvec3 & spherical, cvec3 & cartesian) {
    if (real_t(r)<=EPS){
      std::cout << "Logic Error: encountered spherical to cartesian transformation with r==0!"<<
                      " This is most likely due to L2P with P in center of L"<<"\n";
      throw std::invalid_argument( "Invalid r==0" );
    }
    cartesian[0] = sin(theta) * cos(phi) * spherical[0] // x component (not x itself)
      + cos(theta) * cos(phi) / r * spherical[1]
      - sin(phi) / r / sin(theta) * spherical[2];
    cartesian[1] = sin(theta) * sin(phi) * spherical[0] // y component (not y itself)
      + cos(theta) * sin(phi) / r * spherical[1]
      + cos(phi) / r / sin(theta) * spherical[2];
    cartesian[2] = cos(theta) * spherical[0]                 // z component (not z itself)
      - sin(theta) / r * spherical[1];
  }


  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real_t rho, real_t alpha, real_t beta, complex_t * Ynm, complex_t * YnmTheta) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t invY = y == 0 ? 0 : 1 / y;                           // 1 / y
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t rhom = 1;                                            // Initialize rho^m
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    // NOTE: m and n are switched (Y^n_m instead of Y^m_n)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      // NOTE: The recursion relation is actually
      //                p = ( (2m+1)x*p1 - (m+n)p2 ) / (m-n+1),
      //        so this is only the part needed for the first
      //        derivative. It seems that this is P^{m}_{m+1}
      //        as shown in Eq. 17 in Ref [1].
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim; //  theta derivative of r^n * Ynm
      rhom *= rho;                                              //  rho^m
      real_t rhon = rhom;                                       //  rho^n
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        // NOTE: rhon = (-1)^(m+n_i) * rhom
        //                      / PI^{n_i}_{i=1}(2*m + i),
        //        with n_i the iteration number (n=m+1 => n_i=1)
        //            = (-1)^(n_i) * rho^m
        //                    /  PI^{m}_{i=0}[(2*i+2)*(2*i+1)]
        //                    / PI^{n_i}_{i=1}(2*m + i)
        rhon /= -(n + m);                                       //   Update factorial
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        // NOTE: Here is the recursion applied relation applied
        //        for real. It throws me off that m and n are
        //        not switched anymore by this point.
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        // NOTE: Compare this line with Eq. 23 in Ref. [1]. It
        //        seems that rhon = rho^n * sqrt((n-abs(m))!
        //                                      / (n+abs(m))!)
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;// theta derivative
        rhon *= rho;                                            //   Update rho^n
      }                                                         //  End loop over n in Ynm
      // NOTE: rhom = (-1)^m * rho^m
      //                        /  PI^{m}_{i=0}[(2*i+2)*(2*i+1)]
      rhom /= -(2 * m + 2) * (2 * m + 1);                       //  Update factorial
      // NOTE: pn = (-1)^m * (sin(alpha))^m * PI^{m}_{i=0}(2*i+1)
      pn = -pn * fact * y;                                      //  Pn
      // NOTE: (rhom) * (pn) = { rho^m/PI^{m}_{i=0}
      //                            [(2*i+2)*(2*i+1)] }
      //                  * {(sin(alpha))^m*PI^{m}_{i=0}(2*i+1)}
      //                     = ( rho*sin(alpha) )^m
      //                              / PI^{m}_{i=0}(2*i+2)
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }

  //C^2 multicomplex variable. See notes on 20170915 notebook.
  struct multicomplex {
    cplx A;                                     //!< multicomplex "Real"
    cplx B;                                     //!< multicomplex "Imaginary" (j)
  };
  multicomplex init_multicomplex(cplx A, cplx B){
    multicomplex out;
    out.A = A;
    out.B = B;
    return out;
  }
  multicomplex init_from_C1(complex_t C){
    multicomplex out;
    out.A = cplx(std::real(C), 0);
    out.B = cplx(imag(C), 0);
    return out;
  }
  multicomplex product(multicomplex V, multicomplex W){
      return init_multicomplex(V.A*W.A - V.B*W.B, V.A*W.B + V.B*W.A);
  }
  multicomplex product(multicomplex V, cplx X){
      return init_multicomplex(V.A*X, V.B*X);
  }
  multicomplex product(multicomplex V, real_t x){
      return init_multicomplex(V.A*x, V.B*x);
  }
  // e^{jX} in C^2 where X in C^1
  multicomplex multi_exp(cplx X){
      return init_multicomplex(cos(X), sin(X));
  }
  multicomplex conjugate(multicomplex V){
      return init_multicomplex(V.A, -V.B);
  }
  complex_t Re1(multicomplex V){
    return complex_t(real_t(V.A), real_t(V.B));
  }
  complex_t Im1(multicomplex V){
    return complex_t(imag(V.A), imag(V.B));
  }

  void evalMultipole(cplx rho, cplx alpha, cplx beta, multicomplex * Ynm, multicomplex * YnmTheta) {
    cplx x = cos(alpha);
    cplx y = sin(alpha);
    cplx invY = (real_t(y)<=EPS) ? cplx(0,0) : 1 / y;
    cplx fact = 1;
    cplx pn = 1;
    cplx rhom = 1;
    // complex_t ei = std::exp(I * beta);
    // complex_t eim = 1.0;
    multicomplex ei = multi_exp(beta);
    multicomplex eim = init_multicomplex(1.0, 0.0);
    for (int m=0; m<P; m++) {
      cplx p = pn;
      int npn = m * m + 2 * m;
      int nmn = m * m;
      // Ynm[npn] = rhom * p * eim;
      // Ynm[nmn] = std::conj(Ynm[npn]);
      Ynm[npn] = product(eim, rhom * p);
      Ynm[nmn] = conjugate(Ynm[npn]);
      cplx p1 = p;
      p = x * (2 * m + 1) * p1;
      // YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) * invY * eim;
      YnmTheta[npn] = product(eim, rhom * (p - (m + 1) * x * p1) * invY);
      rhom *= rho;
      cplx rhon = rhom;
      for (int n=m+1; n<P; n++) {
        int npm = n * n + n + m;
        int nmm = n * n + n - m;
        rhon /= -(n + m);
        // Ynm[npm] = rhon * p * eim;
        // Ynm[nmm] = std::conj(Ynm[npm]);
        Ynm[npm] = product(eim, rhon * p);
        Ynm[nmm] = conjugate(Ynm[npm]);
        cplx p2 = p1;
        p1 = p;
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);
        // YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY * eim;
        YnmTheta[npm] = product(eim, rhon * ((n - m + 1) * p - (n + 1) * x * p1) * invY);
        rhon *= rho;
      }
      rhom /= -(2 * m + 2) * (2 * m + 1);
      pn = -pn * fact * y;
      fact += 2;
      // eim *= ei;
      eim = product(eim,ei);
    }
  }


  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real_t rho, real_t alpha, real_t beta, complex_t * Ynm) {
    real_t x = std::cos(alpha);                                 // x = cos(alpha)
    real_t y = std::sin(alpha);                                 // y = sin(alpha)
    real_t fact = 1;                                            // Initialize 2 * m + 1
    real_t pn = 1;                                              // Initialize Legendre polynomial Pn
    real_t invR = -1.0 / rho;                                   // - 1 / rho
    real_t rhom = -invR;                                        // Initialize rho^(-m-1)
    complex_t ei = std::exp(I * beta);                          // exp(i * beta)
    complex_t eim = 1.0;                                        // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                   // Loop over m in Ynm
      real_t p = pn;                                            //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                  //  Index of Ynm for m > 0
      int nmn = m * m;                                          //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                           //  Use conjugate relation for m < 0
      real_t p1 = p;                                            //  Pnm-1
      p = x * (2 * m + 1) * p1;                                 //  Pnm using recurrence relation
      rhom *= invR;                                             //  rho^(-m-1)
      real_t rhon = rhom;                                       //  rho^(-n-1)
      for (int n=m+1; n<P; n++) {                               //  Loop over n in Ynm
        int npm = n * n + n + m;                                //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * eim;                              //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                         //   Use conjugate relation for m < 0
        real_t p2 = p1;                                         //   Pnm-2
        p1 = p;                                                 //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);//   Pnm using recurrence relation
        rhon *= invR * (n - m + 1);                             //   rho^(-n-1)
      }                                                         //  End loop over n in Ynm
      pn = -pn * fact * y;                                      //  Pn
      fact += 2;                                                //  2 * m + 1
      eim *= ei;                                                //  Update exp(i * m * beta)
    }                                                           // End loop over m in Ynm
  }

  void initKernel() {
    NTERM = P * (P + 1) / 2;
  }

  //NOTE: It seems like the P2P implementation flips the sign of the vector
  // potential as seen when it assigns Bi[i].J -= J;, with it doesn't do that
  // when it assigns the potential field p. This required flipping the sign
  // of my derivatives calculations, and I'm not exactly sure how it will turn
  // out when I use the potential field for calculating the Gaussian of the RBF.

  // Gaussian-regularized kernel that calculates the Hessian analytically
  void P2P_reggausANA(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      vec9 J = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR2 = 1.0 / R2;
          real_t invR3 = 1.0 / R2 / sqrt(R2);
          real_t sgm2 = Bj[j].sigma[0]*Bj[j].sigma[0];
          real_t xi = R2/2/sgm2;
          real_t sqrtxi = sqrt(xi);
          real_t gsgm = erf(sqrtxi) - SQRT4pi*sqrtxi*std::exp(-xi);
          vec3 dgsgmdx = dX*SQRT4pi*sqrtxi*std::exp(-xi)/sgm2;
          vec3 K = dX*invR3;                            // Singular kernel
          vec3 Ksgm = K*gsgm;                           // Regularized kernel
          for(int ind=0; ind<3; ind++){
            vec3 aux1 = Ksgm*Bj[j].q[ind];
            vec3 aux2 = dgsgmdx*Bj[j].q[ind];
            real_t aux3 = -3*gsgm*Bj[j].q[ind]*invR3*invR2;
            // p[ind] += qinvR;                         // No need for calculating the potential
            // Jacobian J[i*3 + j] = dFi/dxj
            J[ind*3 + 0] += aux1[0];
            J[ind*3 + 1] += aux1[1];
            J[ind*3 + 2] += aux1[2];
            // Hessian dJdxk[i*3 + j] = d(dFi/dxj)/dxk
            Bi[i].dJdx1[ind*3 + 0] -= aux2[0]*K[0] + aux3*dX[0]*dX[0] + gsgm*Bj[j].q[ind]*invR3;
            Bi[i].dJdx1[ind*3 + 1] -= aux2[0]*K[1] + aux3*dX[1]*dX[0];
            Bi[i].dJdx1[ind*3 + 2] -= aux2[0]*K[2] + aux3*dX[2]*dX[0];
            Bi[i].dJdx2[ind*3 + 0] -= aux2[1]*K[0] + aux3*dX[0]*dX[1];
            Bi[i].dJdx2[ind*3 + 1] -= aux2[1]*K[1] + aux3*dX[1]*dX[1] + gsgm*Bj[j].q[ind]*invR3;
            Bi[i].dJdx2[ind*3 + 2] -= aux2[1]*K[2] + aux3*dX[2]*dX[1];
            Bi[i].dJdx3[ind*3 + 0] -= aux2[2]*K[0] + aux3*dX[0]*dX[2];
            Bi[i].dJdx3[ind*3 + 1] -= aux2[2]*K[1] + aux3*dX[1]*dX[2];
            Bi[i].dJdx3[ind*3 + 2] -= aux2[2]*K[2] + aux3*dX[2]*dX[2] + gsgm*Bj[j].q[ind]*invR3;
          }
        }
      }
      Bi[i].p += p;
      Bi[i].J -= J;
    }
  }



  // Singular kernel
  void P2P_sing(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      vec9 J = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR2 = 1.0 / R2;
          for(int ind=0; ind<3; ind++){
            real_t invR = Bj[j].q[ind] * sqrt(invR2);
            p[ind] += invR;
            vec3 aux1 = dX * invR2 * invR;
            J[ind*3 + 0] += aux1[0];
            J[ind*3 + 1] += aux1[1];
            J[ind*3 + 2] += aux1[2];
          }
        }
      }
      Bi[i].p += p;
      Bi[i].J -= J;
    }
  }


  // Singular kernel that calculates the Hessian analytically
  void P2P_singANA(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      vec9 J = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t invR2 = 1.0 / R2;
          real_t invR = sqrt(invR2);
          real_t invR3 = invR*invR2;
          real_t invR4 = invR2*invR2;
          for(int ind=0; ind<3; ind++){
            real_t qinvR = Bj[j].q[ind] * invR;
            vec3 aux1 = dX * invR2 * qinvR;
            // Potential p[i] = Fi
            p[ind] += qinvR;
            // Jacobian J[i*3 + j] = dFi/dxj
            J[ind*3 + 0] += aux1[0];
            J[ind*3 + 1] += aux1[1];
            J[ind*3 + 2] += aux1[2];
            // Hessian dJdxk[i*3 + j] = d(dFi/dxj)/dxk
            Bi[i].dJdx1[ind*3 + 0] += 3*qinvR*invR4*dX[0]*dX[0] - Bj[j].q[ind]*invR3;
            Bi[i].dJdx1[ind*3 + 1] += 3*qinvR*invR4*dX[1]*dX[0];
            Bi[i].dJdx1[ind*3 + 2] += 3*qinvR*invR4*dX[2]*dX[0];
            Bi[i].dJdx2[ind*3 + 0] += 3*qinvR*invR4*dX[0]*dX[1];
            Bi[i].dJdx2[ind*3 + 1] += 3*qinvR*invR4*dX[1]*dX[1] - Bj[j].q[ind]*invR3;
            Bi[i].dJdx2[ind*3 + 2] += 3*qinvR*invR4*dX[2]*dX[1];
            Bi[i].dJdx3[ind*3 + 0] += 3*qinvR*invR4*dX[0]*dX[2];
            Bi[i].dJdx3[ind*3 + 1] += 3*qinvR*invR4*dX[1]*dX[2];
            Bi[i].dJdx3[ind*3 + 2] += 3*qinvR*invR4*dX[2]*dX[2] - Bj[j].q[ind]*invR3;
          }
        }
      }
      Bi[i].p += p;
      Bi[i].J -= J;
    }
  }

  // Regularized (Winckelmann's) kernel
  void P2P_reg(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      vec9 J = 0;
      vec3 pse = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t sgm2 = Bj[j].sigma[0]*Bj[j].sigma[0];
          real_t aux2 = R2 / (sgm2);
          real_t aux3 = (aux2 + 1.5) / std::pow(aux2 + 1, 1.5) / Bj[j].sigma[0];
          real_t aux4 = (aux2 + 2.5) / std::pow(aux2 + 1, 2.5) / (sgm2*Bj[j].sigma[0]);
          vec3 aux1 = dX * aux4;
          for(int ind=0; ind<3; ind++){
            p[ind] += Bj[j].q[ind] * aux3;
            J[ind*3 + 0] += aux1[0] * Bj[j].q[ind];
            J[ind*3 + 1] += aux1[1] * Bj[j].q[ind];
            J[ind*3 + 2] += aux1[2] * Bj[j].q[ind];
          }
          // Particle Strength Exchange
          real_t sgmij2 = (Bi[i].sigma[0]*Bi[i].sigma[0] + Bj[j].sigma[0]*Bj[j].sigma[0])/2;
          real_t dnmntr = std::pow(sgmij2, 2.5) * std::pow(R2/sgmij2 + 1, 4.5);
          for(int ind=0; ind<3; ind++){
            pse[ind] += (Bi[i].vol[0]*Bj[j].q[ind] - Bj[j].vol[0]*Bi[i].q[ind])/dnmntr;
          }
        }
      }
      Bi[i].p += p;
      Bi[i].J -= J;
      Bi[i].pse += pse;
    }
  }

  // Regularized (Winckelmann's) kernel that calculates the Hessian analytically
  void P2P_regANA(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;
      vec9 J = 0;
      vec3 pse = 0;
      for (int j=0; j<Cj->numBodies; j++) {
        vec3 dX = Bi[i].X - Bj[j].X;
        real_t R2 = norm(dX);
        if (R2 != 0) {
          real_t sgm2 = Bj[j].sigma[0]*Bj[j].sigma[0];
          real_t aux2 = R2 / (sgm2);
          real_t aux3 = (aux2 + 1.5) / std::pow(aux2 + 1, 1.5) / Bj[j].sigma[0];
          real_t aux4 = (aux2 + 2.5) / std::pow(aux2 + 1, 2.5) / (sgm2*Bj[j].sigma[0]);
          real_t aux5 = (3*aux2+10.5) / std::pow(aux2 + 1, 3.5) / (sgm2*sgm2*Bj[j].sigma[0]);
          vec3 aux1 = dX * aux4;

          // std::ostringstream oss;
          // oss << "Jexa_ncrit" << NCRIT << "-" << Bi[i].index[0] << "-" << Bj[j].index[0] << ".csv";
          // std::string fname = oss.str();
          // ofstream myfile;
          // myfile.open(fname);
          // myfile.precision(12);
          // for(int ind=0; ind<3; ind++){
          //   myfile << aux1[0] * Bj[j].q[ind] << " ";
          //   myfile << aux1[1] * Bj[j].q[ind] << " ";
          //   myfile << aux1[2] * Bj[j].q[ind];
          //   if(ind!=2) myfile << " ";
          // }
          // myfile.close();

          for(int ind=0; ind<3; ind++){
          // Potential p[i] = Fi
            p[ind] += Bj[j].q[ind] * aux3;
            // Jacobian J[i*3 + j] = dFi/dxj
            J[ind*3 + 0] += aux1[0] * Bj[j].q[ind];
            J[ind*3 + 1] += aux1[1] * Bj[j].q[ind];
            J[ind*3 + 2] += aux1[2] * Bj[j].q[ind];
            // Hessian dJdxk[i*3 + j] = d(dFi/dxj)/dxk
            Bi[i].dJdx1[ind*3 + 0] += Bj[j].q[ind]*(dX[0]*dX[0]*aux5 - aux4);
            Bi[i].dJdx1[ind*3 + 1] += Bj[j].q[ind]*(dX[1]*dX[0]*aux5);
            Bi[i].dJdx1[ind*3 + 2] += Bj[j].q[ind]*(dX[2]*dX[0]*aux5);
            Bi[i].dJdx2[ind*3 + 0] += Bj[j].q[ind]*(dX[0]*dX[1]*aux5);
            Bi[i].dJdx2[ind*3 + 1] += Bj[j].q[ind]*(dX[1]*dX[1]*aux5 - aux4);
            Bi[i].dJdx2[ind*3 + 2] += Bj[j].q[ind]*(dX[2]*dX[1]*aux5);
            Bi[i].dJdx3[ind*3 + 0] += Bj[j].q[ind]*(dX[0]*dX[2]*aux5);
            Bi[i].dJdx3[ind*3 + 1] += Bj[j].q[ind]*(dX[1]*dX[2]*aux5);
            Bi[i].dJdx3[ind*3 + 2] += Bj[j].q[ind]*(dX[2]*dX[2]*aux5 - aux4);
          }
          real_t sgmij2 = (Bi[i].sigma[0]*Bi[i].sigma[0] + Bj[j].sigma[0]*Bj[j].sigma[0])/2;
          real_t dnmntr = std::pow(sgmij2, 2.5) * std::pow(R2/sgmij2 + 1, 4.5);
          for(int ind=0; ind<3; ind++){
            pse[ind] += (Bi[i].vol[0]*Bj[j].q[ind] - Bj[j].vol[0]*Bi[i].q[ind])/dnmntr;
          }
        }
      }
      Bi[i].p += p;
      Bi[i].J -= J;
      Bi[i].pse += pse;
    }
  }



  // Regularized (Winckelmann's) kernel that calculates the Hessian through CSDA
  void P2P_regCSDA(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      for(int H_ind=0; H_ind<3; H_ind++){
        cvec3 p;
        cvec9 J;
        for(int ind=0; ind<3; ind++) p[ind] = 0;
        for(int ind=0; ind<9; ind++) J[ind] = 0;
        for (int j=0; j<Cj->numBodies; j++) {
          cvec3 dX;
          for(int ind=0; ind<3; ind++) dX[ind] = Bi[i].X[ind] - Bj[j].X[ind];
          dX[H_ind] += cplx(0, COMPLEX_STEP);
          cplx R2 = norm(dX);
          if (abs(real_t(R2)) > EPS) {
            real_t sgm2 = Bj[j].sigma[0]*Bj[j].sigma[0];
            cplx aux2 = R2 / (sgm2);
            cplx aux3 = (aux2 + real_t(1.5)) / pow(aux2 + 1, real_t(1.5)) / Bj[j].sigma[0];
            cplx aux4 = (aux2 + real_t(2.5)) / pow(aux2 + 1, real_t(2.5)) / (sgm2*Bj[j].sigma[0]);
            cvec3 aux1;
            for(int ind=0; ind<3; ind++) aux1[ind] = dX[ind] * aux4;
            for(int ind=0; ind<3; ind++){
              p[ind] += Bj[j].q[ind] * aux3;
              J[ind*3 + 0] += aux1[0] * Bj[j].q[ind];
              J[ind*3 + 1] += aux1[1] * Bj[j].q[ind];
              J[ind*3 + 2] += aux1[2] * Bj[j].q[ind];
            }
          }
        }
        if(H_ind==0){
          for(int ind=0; ind<3; ind++) Bi[i].p[ind] += real_t(p[ind]);
          for(int ind=0; ind<9; ind++) Bi[i].J[ind] -= real_t(J[ind]);
          for(int ind=0; ind<9; ind++) Bi[i].dJdx1[ind] -= imag(J[ind])/COMPLEX_STEP;
        }
        else if(H_ind==1){
          for(int ind=0; ind<9; ind++) Bi[i].dJdx2[ind] -= imag(J[ind])/COMPLEX_STEP;
        }
        else if(H_ind==2){
          for(int ind=0; ind<9; ind++) Bi[i].dJdx3[ind] -= imag(J[ind])/COMPLEX_STEP;
        }
      }
    }
  }

  void P2M_a(Cell * C) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=C->body; B!=C->body+C->numBodies; B++) {
      vec3 dX = B->X - C->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, -beta, Ynm, YnmTheta);
      for (int n=0; n<P; n++) {
        for (int m=0; m<=n; m++) {
          int nm  = n * n + n + m;
          int nms = n * (n + 1) / 2 + m;
          for (int i=0; i<3; i++) C->M[i][nms] += B->q[i] * Ynm[nm];
        }
      }
    }
  }

  void M2M_a(Cell * Ci) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Cell * Cj=Ci->child; Cj!=Ci->child+Ci->numChilds; Cj++) {
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          std::vector<complex_t> M;
          M.resize(3, 0.0);
          for (int n=0; n<=j; n++) {
            for (int m=std::max(-n,-j+k+n); m<=std::min(k-1,n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
              int nm    = n * n + n - m;
              for(int ind=0; ind<3; ind++){
                M[ind] += Cj->M[ind][jnkms] * Ynm[nm] * real_t(ipow2n(m) * oddOrEven(n));
              }
            }
            for (int m=k; m<=std::min(n,j+k-n); m++) {
              int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
              int nm    = n * n + n - m;
              for(int ind=0; ind<3; ind++){
                M[ind] += std::conj(Cj->M[ind][jnkms]) * Ynm[nm] * real_t(oddOrEven(k+n+m));
              }
            }
          }
          for(int ind=0; ind<3; ind++) Ci->M[ind][jks] += M[ind];
        }
      }
    }
  }

  void M2L_a(Cell * Ci, Cell * Cj) {
    complex_t Ynm2[4*P*P];
    vec3 dX = Ci->X - Cj->X;
    real_t rho, alpha, beta;
    cart2sph(dX, rho, alpha, beta);
    evalLocal(rho, alpha, beta, Ynm2);
    for (int j=0; j<P; j++) {
      real_t Cnm = oddOrEven(j);
      for (int k=0; k<=j; k++) {
        int jks = j * (j + 1) / 2 + k;
        std::vector<complex_t> L;
        L.resize(3, 0.0);
        for (int n=0; n<P; n++) {
          for (int m=-n; m<0; m++) {
            int nms  = n * (n + 1) / 2 - m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            for(int ind=0; ind<3; ind++) L[ind] += std::conj(Cj->M[ind][nms]) * Cnm * Ynm2[jnkm];
          }
          for (int m=0; m<=n; m++) {
            int nms  = n * (n + 1) / 2 + m;
            int jnkm = (j + n) * (j + n) + j + n + m - k;
            real_t Cnm2 = Cnm * oddOrEven((k-m)*(k<m)+m);
            for(int ind=0; ind<3; ind++) L[ind] += Cj->M[ind][nms] * Cnm2 * Ynm2[jnkm];
          }
        }
        for(int ind=0; ind<3; ind++) Ci->L[ind][jks] += L[ind];
      }
    }
  }

  void L2L_a(Cell * Cj) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
      vec3 dX = Ci->X - Cj->X;
      real_t rho, alpha, beta;
      cart2sph(dX, rho, alpha, beta);
      evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
      for (int j=0; j<P; j++) {
        for (int k=0; k<=j; k++) {
          int jks = j * (j + 1) / 2 + k;
          std::vector<complex_t> L;
          L.resize(3, 0.0);
          for (int n=j; n<P; n++) {
            for (int m=j+k-n; m<0; m++) {
              int jnkm = (n - j) * (n - j) + n - j + m - k;
              int nms  = n * (n + 1) / 2 - m;
              for(int ind=0; ind<3; ind++){
                L[ind] += std::conj(Cj->L[ind][nms]) * Ynm[jnkm] * real_t(oddOrEven(k));
              }
            }
            for (int m=0; m<=n; m++) {
              if (n-j >= abs(m-k)) {
                int jnkm = (n - j) * (n - j) + n - j + m - k;
                int nms  = n * (n + 1) / 2 + m;
                for(int ind=0; ind<3; ind++){
                  L[ind] += Cj->L[ind][nms] * Ynm[jnkm] * real_t(oddOrEven((m-k)*(m<k)));
                }
              }
            }
          }
          for(int ind=0; ind<3; ind++) Ci->L[ind][jks] += L[ind];
        }
      }
    }
  }

  // L2P without calculating the Hessian
  void L2P_noH(Cell * Ci) {
    complex_t Ynm[P*P], YnmTheta[P*P];
    for (Body * B=Ci->body; B!=Ci->body+Ci->numBodies; B++) {
      vec3 dX = B->X - Ci->X;
      vec3 spherical1 = 0;
      vec3 spherical2 = 0;
      vec3 spherical3 = 0;
      vec3 cartesian = 0;
      real_t r, theta, phi;
      cart2sph(dX, r, theta, phi);
      evalMultipole(r, theta, phi, Ynm, YnmTheta);
      for (int n=0; n<P; n++) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        for(int ind=0; ind<3; ind++) B->p[ind] += std::real(Ci->L[ind][nms] * Ynm[nm]);
        spherical1[0] += std::real(Ci->L[0][nms] * Ynm[nm]) / r * n;
        spherical1[1] += std::real(Ci->L[0][nms] * YnmTheta[nm]);
        spherical2[0] += std::real(Ci->L[1][nms] * Ynm[nm]) / r * n;
        spherical2[1] += std::real(Ci->L[1][nms] * YnmTheta[nm]);
        spherical3[0] += std::real(Ci->L[2][nms] * Ynm[nm]) / r * n;
        spherical3[1] += std::real(Ci->L[2][nms] * YnmTheta[nm]);
        for (int m=1; m<=n; m++) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          for(int ind=0; ind<3; ind++) B->p[ind] += 2 * std::real(Ci->L[ind][nms] * Ynm[nm]);
          spherical1[0] += 2 * std::real(Ci->L[0][nms] * Ynm[nm]) / r * n;
          spherical1[1] += 2 * std::real(Ci->L[0][nms] * YnmTheta[nm]);
          spherical1[2] += 2 * std::real(Ci->L[0][nms] * Ynm[nm] * I) * m;
          spherical2[0] += 2 * std::real(Ci->L[1][nms] * Ynm[nm]) / r * n;
          spherical2[1] += 2 * std::real(Ci->L[1][nms] * YnmTheta[nm]);
          spherical2[2] += 2 * std::real(Ci->L[1][nms] * Ynm[nm] * I) * m;
          spherical3[0] += 2 * std::real(Ci->L[2][nms] * Ynm[nm]) / r * n;
          spherical3[1] += 2 * std::real(Ci->L[2][nms] * YnmTheta[nm]);
          spherical3[2] += 2 * std::real(Ci->L[2][nms] * Ynm[nm] * I) * m;
        }
      }
      cartesian = 0;
      sph2cart(r, theta, phi, spherical1, cartesian);
      for(int ind=0; ind<3; ind++) B->J[3*0 + ind] += cartesian[ind];
      cartesian = 0;
      sph2cart(r, theta, phi, spherical2, cartesian);
      for(int ind=0; ind<3; ind++) B->J[3*1 + ind] += cartesian[ind];
      cartesian = 0;
      sph2cart(r, theta, phi, spherical3, cartesian);
      for(int ind=0; ind<3; ind++) B->J[3*2 + ind] += cartesian[ind];
    }
  }

  // L2P calculating the Hessian through CSDA
  void L2P_CSDA(Cell * Ci) {
    multicomplex Ynm[P*P], YnmTheta[P*P];
    for (Body * B=Ci->body; B!=Ci->body+Ci->numBodies; B++) {
      for(int H_ind=0; H_ind<3; H_ind++){
        // vec3 dX = B->X - Ci->X;
        cvec3 dX;
        for(int ind=0; ind<3; ind++) dX[ind] = B->X[ind] - Ci->X[ind];
        dX[H_ind] += cplx(0, COMPLEX_STEP);
        cvec3 spherical1;
        cvec3 spherical2;
        cvec3 spherical3;
        cvec3 cartesian;
        for(int ind=0; ind<3; ind++){
          spherical1[ind] = 0;
          spherical2[ind] = 0;
          spherical3[ind] = 0;
          cartesian[ind] = 0;
        }
        cplx r, theta, phi;
        cart2sph(dX, r, theta, phi);
        evalMultipole(r, theta, phi, Ynm, YnmTheta);
        for (int n=0; n<P; n++) {
          int nm  = n * n + n;
          int nms = n * (n + 1) / 2;
          if(H_ind==0){
            for(int ind=0; ind<3; ind++){
              B->p[ind] += real_t((product(init_from_C1(Ci->L[ind][nms]), Ynm[nm])).A);}
          }
          spherical1[0] += (product(init_from_C1(Ci->L[0][nms]), Ynm[nm])).A / r * n;
          spherical1[1] += (product(init_from_C1(Ci->L[0][nms]), YnmTheta[nm])).A;
          spherical2[0] += (product(init_from_C1(Ci->L[1][nms]), Ynm[nm])).A / r * n;
          spherical2[1] += (product(init_from_C1(Ci->L[1][nms]), YnmTheta[nm])).A;
          spherical3[0] += (product(init_from_C1(Ci->L[2][nms]), Ynm[nm])).A / r * n;
          spherical3[1] += (product(init_from_C1(Ci->L[2][nms]), YnmTheta[nm])).A;
          for (int m=1; m<=n; m++) {
            nm  = n * n + n + m;
            nms = n * (n + 1) / 2 + m;
            if(H_ind==0){
              for(int ind=0; ind<3; ind++){
                B->p[ind] += real_t(2 * (product(init_from_C1(Ci->L[ind][nms]), Ynm[nm])).A);}
            }
            spherical1[0] += 2 * (product(init_from_C1(Ci->L[0][nms]), Ynm[nm])).A / r * n;
            spherical1[1] += 2 * (product(init_from_C1(Ci->L[0][nms]), YnmTheta[nm])).A;
            spherical1[2] += 2 * (product( product(init_from_C1(Ci->L[0][nms]), Ynm[nm]), init_from_C1(I) )).A * m;
            spherical2[0] += 2 * (product(init_from_C1(Ci->L[1][nms]), Ynm[nm])).A / r * n;
            spherical2[1] += 2 * (product(init_from_C1(Ci->L[1][nms]), YnmTheta[nm])).A;
            spherical2[2] += 2 * (product( product(init_from_C1(Ci->L[1][nms]), Ynm[nm]), init_from_C1(I) )).A * m;
            spherical3[0] += 2 * (product(init_from_C1(Ci->L[2][nms]), Ynm[nm])).A / r * n;
            spherical3[1] += 2 * (product(init_from_C1(Ci->L[2][nms]), YnmTheta[nm])).A;
            spherical3[2] += 2 * (product( product(init_from_C1(Ci->L[2][nms]), Ynm[nm]), init_from_C1(I) )).A * m;
          }
        }
        for(int ind=0; ind<3; ind++) cartesian[ind] = 0;
        sph2cart(r, theta, phi, spherical1, cartesian);

        if(H_ind==0){
          for(int ind=0; ind<3; ind++){
            B->J[3*0 + ind] += real_t(cartesian[ind]);
            B->dJdx1[3*0 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==1){
          for(int ind=0; ind<3; ind++){
            B->dJdx2[3*0 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==2){
          for(int ind=0; ind<3; ind++){
            B->dJdx3[3*0 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }

        for(int ind=0; ind<3; ind++) cartesian[ind] = 0;
        sph2cart(r, theta, phi, spherical2, cartesian);

        if(H_ind==0){
          for(int ind=0; ind<3; ind++){
            B->J[3*1 + ind] += real_t(cartesian[ind]);
            B->dJdx1[3*1 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==1){
          for(int ind=0; ind<3; ind++){
            B->dJdx2[3*1 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==2){
          for(int ind=0; ind<3; ind++){
            B->dJdx3[3*1 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }

        for(int ind=0; ind<3; ind++) cartesian[ind] = 0;
        sph2cart(r, theta, phi, spherical3, cartesian);

        if(H_ind==0){
          for(int ind=0; ind<3; ind++){
            B->J[3*2 + ind] += real_t(cartesian[ind]);
            B->dJdx1[3*2 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==1){
          for(int ind=0; ind<3; ind++){
            B->dJdx2[3*2 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
        else if(H_ind==2){
          for(int ind=0; ind<3; ind++){
            B->dJdx3[3*2 + ind] += imag(cartesian[ind])/COMPLEX_STEP;
          }
        }
      }

    }
  }


  // Gaussian field discretization kernel
  void GaussianzetaP2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;

      for (int j=0; j<Cj->numBodies; j++) {
        real_t R2 = norm(Bi[i].X - Bj[j].X);

        // if (R2 != 0) {
          real_t S2 = 2 * Bj[j].sigma[0] * Bj[j].sigma[0];
          real_t aux1 = std::exp(-R2 / S2) / (M_PI * S2) / sqrt(M_PI * S2);
            p += Bj[j].q * aux1;
        // }

      }

      Bi[i].J[0] += p[0];
      Bi[i].J[1] += p[1];
      Bi[i].J[2] += p[2];
    }
  }


  // Algebraic field discretization kernel
  void algzetaP2P(Cell * Ci, Cell * Cj) {
    Body * Bi = Ci->body;
    Body * Bj = Cj->body;
    for (int i=0; i<Ci->numBodies; i++) {
      vec3 p = 0;

      for (int j=0; j<Cj->numBodies; j++) {
        real_t R2 = norm(Bi[i].X - Bj[j].X);

        // if (R2 != 0) {
          real_t aux1 = R2 / (Bj[j].sigma[0]*Bj[j].sigma[0]);
          real_t aux2 = 7.5 / std::pow(aux1 + 1, 3.5) / (Bj[j].sigma[0]*Bj[j].sigma[0]*Bj[j].sigma[0]);
          p += Bj[j].q * aux2 / (M_PI*4);
        // }

      }

      Bi[i].J[0] += p[0];
      Bi[i].J[1] += p[1];
      Bi[i].J[2] += p[2];
    }
  }


  void P2P_a(Cell * Ci, Cell * Cj){
    switch(P2P_TYPE){
      case 0: P2P_sing(Ci, Cj); break;
      case 1: P2P_singANA(Ci, Cj); break;
      case 2: P2P_reg(Ci, Cj); break;
      case 3: P2P_regANA(Ci, Cj); break;
      case 4: P2P_regCSDA(Ci, Cj); break;
      case 5: P2P_reggausANA(Ci, Cj); break;
      default: std::cout << "Exception 123: Invalid P2P_TYPE "<<P2P_TYPE<<"\n"; throw 123;
    }
  }
  void L2P_a(Cell * Ci){
    switch(L2P_TYPE){
      case 0: L2P_noH(Ci); break;
      case 1: L2P_CSDA(Ci); break;
      default: std::cout << "Exception 124: Invalid L2P_TYPE "<<L2P_TYPE<<"\n"; throw 124;
    }
  }

  void zetaP2P(Cell * Ci, Cell * Cj){
    switch(P2P_TYPE){
      case 3: algzetaP2P(Ci, Cj); break;
      case 5: GaussianzetaP2P(Ci, Cj); break;
      default: std::cout << "Exception 125: Invalid P2P_TYPE "<<P2P_TYPE<<"\n"; throw 123;
    }
  }


  void P2P(Cell * Ci, Cell * Cj){ if(RBF){ zetaP2P(Ci, Cj); }
                                  else{ P2P_a(Ci, Cj); } }
  void L2P(Cell * Ci){ if(RBF==false) L2P_a(Ci); }
  void P2M(Cell * C){ if(RBF==false) P2M_a(C); }
  void M2M(Cell * Ci){ if(RBF==false) M2M_a(Ci); }
  void M2L(Cell * Ci, Cell * Cj){ if(RBF==false) M2L_a(Ci, Cj); }
  void L2L(Cell * Cj){ if(RBF==false) L2L_a(Cj); }


}
#endif
