#ifndef exafmm_h
#define exafmm_h
#include <complex>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cmath>
#include "vec.h"
#include "complexify.h"


namespace exafmm {
  //! Basic type definitions
#if EXAFMM_SINGLE
  typedef float real_t;                         //!< Floating point type is single precision
  const real_t EPS = 1e-8f;                     //!< Single precision epsilon
  const real_t COMPLEX_STEP = 1e-8f;            //!< Complex step for derivative approximation
  const bool single_prec = true;
#else
  typedef double real_t;                        //!< Floating point type is double precision
  const real_t EPS = 1e-16;                     //!< Double precision epsilon
  const real_t COMPLEX_STEP = 1e-16;            //!< Complex step for derivative approximation
  const bool single_prec = false;
#endif
  typedef std::complex<real_t> complex_t;       //!< Complex type
  typedef vec<1,real_t> vec1;                   //!< Vector of 1 real_t types
  typedef vec<1,int> vec1int;                   //!< Vector of 1 int types
  typedef vec<3,real_t> vec3;                   //!< Vector of 3 real_t types
  typedef vec<3,cplx> cvec3;                    //!< Vector of 3 cplx types
  typedef vec<9,real_t> vec9;                   //!< Vector of 9 real_t types
  typedef vec<9,cplx> cvec9;                    //!< Vector of 9 cplx types
  typedef vec<2,cplx> cplxofcplx;               //!< The complex number of a complex

  //! Structure of bodies
  struct Body {
    // Properties
    vec3 X;                                     //!< Position
    vec3 q;                                     //!< Vector charge
    vec3 vort;                                  //!< Vorticity for RBF
    vec1 sigma;                                 //!< Vortex blob radius
    vec1 vol;                                   //!< Vortex particle volume (it will also use this for RBF)
    vec1int index;                              //!< Index of this body
    // Calculations
    vec3 p;                                     //!< Vector potential
    vec9 J;                                     //!< Jacobian
    vec9 dJdx1;                                 //!< Partial of Jacobian to x1
    vec9 dJdx2;                                 //!< Partial of Jacobian to x2
    vec9 dJdx3;                                 //!< Partial of Jacobian to x3
    vec3 pse;                                   //!< Particle Strength Exchange
    // real_t q;                                   //!< Charge
    // real_t p;                                   //!< Potential
    // vec3 F;                                     //!< Force

    // Functions for exposing internal variables to Julia (allowing to modify them from Julia)
    real_t* get_Xref_aux() { return this->X; }
    jlcxx::ArrayRef<real_t, 1> get_Xref() { return jlcxx::ArrayRef<real_t>(get_Xref_aux(), 3); }
    real_t get_Xi(const int i) { return this->X[i]; }
    void set_Xi(const int i, const real_t val) { this->X[i] = val; }

    real_t* get_qref_aux() { return this->q; }
    jlcxx::ArrayRef<real_t, 1> get_qref() { return jlcxx::ArrayRef<real_t>(get_qref_aux(), 3); }
    real_t get_qi(const int i) { return this->q[i]; }
    void set_qi(const int i, const real_t val) { this->q[i] = val; }

    real_t* get_vortref_aux() { return this->vort; }
    jlcxx::ArrayRef<real_t, 1> get_vortref() { return jlcxx::ArrayRef<real_t>(get_vortref_aux(), 3); }
    real_t get_vorti(const int i) { return this->vort[i]; }
    void set_vorti(const int i, const real_t val) { this->vort[i] = val; }

    real_t* get_sigmaref_aux() { return this->sigma; }
    jlcxx::ArrayRef<real_t, 1> get_sigmaref() { return jlcxx::ArrayRef<real_t>(get_sigmaref_aux(), 1); }
    real_t get_sigmai(const int i) { return this->sigma[i]; }
    void set_sigmai(const int i, const real_t val) { this->sigma[i] = val; }

    real_t* get_volref_aux() { return this->vol; }
    jlcxx::ArrayRef<real_t, 1> get_volref() { return jlcxx::ArrayRef<real_t>(get_volref_aux(), 1); }
    real_t get_voli(const int i) { return this->vol[i]; }
    void set_voli(const int i, const real_t val) { this->vol[i] = val; }

    int* get_indexref_aux() { return this->index; }
    jlcxx::ArrayRef<int, 1> get_indexref() { return jlcxx::ArrayRef<int>(get_indexref_aux(), 1); }
    int get_indexi(const int i) { return this->index[i]; }
    void set_indexi(const int i, const int val) { this->index[i] = val; }

    real_t* get_pref_aux() { return this->p; }
    jlcxx::ArrayRef<real_t, 1> get_pref() { return jlcxx::ArrayRef<real_t>(get_pref_aux(), 3); }
    real_t get_pi(const int i) { return this->p[i]; }
    void set_pi(const int i, const real_t val) { this->p[i] = val; }

    real_t* get_Jref_aux() { return this->J; }
    jlcxx::ArrayRef<real_t, 2> get_Jref() { return jlcxx::ArrayRef<real_t, 2>(get_Jref_aux(), 3, 3); }
    real_t get_Ji(const int i) { return this->J[i]; }
    void set_Ji(const int i, const real_t val) { this->J[i] = val; }

    real_t* get_dJdx1ref_aux() { return this->dJdx1; }
    jlcxx::ArrayRef<real_t, 2> get_dJdx1ref() { return jlcxx::ArrayRef<real_t, 2>(get_dJdx1ref_aux(), 3, 3); }
    real_t get_dJdx1i(const int i) { return this->dJdx1[i]; }
    void set_dJdx1i(const int i, const real_t val) { this->dJdx1[i] = val; }

    real_t* get_dJdx2ref_aux() { return this->dJdx2; }
    jlcxx::ArrayRef<real_t, 2> get_dJdx2ref() { return jlcxx::ArrayRef<real_t, 2>(get_dJdx2ref_aux(), 3, 3); }
    real_t get_dJdx2i(const int i) { return this->dJdx2[i]; }
    void set_dJdx2i(const int i, const real_t val) { this->dJdx2[i] = val; }

    real_t* get_dJdx3ref_aux() { return this->dJdx3; }
    jlcxx::ArrayRef<real_t, 2> get_dJdx3ref() { return jlcxx::ArrayRef<real_t, 2>(get_dJdx3ref_aux(), 3, 3); }
    real_t get_dJdx3i(const int i) { return this->dJdx3[i]; }
    void set_dJdx3i(const int i, const real_t val) { this->dJdx3[i] = val; }

    real_t* get_pseref_aux() { return this->pse; }
    jlcxx::ArrayRef<real_t, 1> get_pseref() { return jlcxx::ArrayRef<real_t>(get_pseref_aux(), 3); }
    real_t get_psei(const int i) { return this->pse[i]; }
    void set_psei(const int i, const real_t val) { this->pse[i] = val; }
  };
  typedef std::vector<Body> Bodies;             //!< Vector of bodies

  //! Structure of cells
  struct Cell {
    int numChilds;                              //!< Number of child cells
    int numBodies;                              //!< Number of descendant bodies
    Cell * child;                               //!< Pointer of first child cell
    Body * body;                                //!< Pointer of first body
    vec3 X;                                     //!< Cell center
    real_t R;                                   //!< Cell radius
    real_t sigma;                               //!< Cell characteristic smoothing radius
#if EXAFMM_LAZY
    std::vector<Cell*> listM2L;                 //!< M2L interaction list
    std::vector<Cell*> listP2P;                 //!< P2P interaction list
#endif
    // std::vector<complex_t> M;                   //!< Multipole expansion coefs
    // std::vector<complex_t> L;                   //!< Local expansion coefs
    std::vector< std::vector<complex_t> > M;    //!< Multipole expansion vectorized coefs
    std::vector< std::vector<complex_t> > L;    //!< Local expansion vectorized coefs
  };
  typedef std::vector<Cell> Cells;              //!< Vector of cells

  //! Global variables
  int P;                                        //!< Order of expansions
  int NTERM;                                    //!< Number of coefficients
  int NCRIT;                                    //!< Number of bodies per leaf cell
  real_t THETA;                                 //!< Multipole acceptance criterion
  real_t PHI;                                   //!< Multipole regularizing acceptance criterion
  int P2P_TYPE=0;                               //!< Flag for the P2P function to call
  int L2P_TYPE=0;                               //!< Flag for the L2P function to call
  // bool CSDA;                                    //!< Complex-step derivative approximation
  // bool REGULARIZED;                             //!< Regularized kernel (and no CSDA)
  real_t SQRT4pi = sqrt(1/M_PI_4);
  real_t LAMBDA = 1.5;                          //! Target core overlap when fixing Lagrangian distortion
  bool RBF = false;                             //! Flag for running RBF routine
  bool SGS = false;                             //! Flag for running SGS routine
  int SGS_TYPE = -1;                            //! Flag for running SGS routine
  bool TRANSPOSED = true;                       //! Flag for using transposed SGS scheme
}
#endif
