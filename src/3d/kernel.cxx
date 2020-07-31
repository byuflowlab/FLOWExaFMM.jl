#include "args.h"
#include "kernel.h"
#include "verify.h"
using namespace exafmm;

int main(int argc, char ** argv) {
  Args args(argc, argv);
  P = args.P;
  // std::cout << "Rabbit2 P=" << P << "\targs.P=" << args.P << "\n";
  VERBOSE = args.verbose;
  initKernel();

  //! P2M
  Bodies jbodies(1);
  jbodies[0].X = 2;
  jbodies[0].q[0] = 1;
  jbodies[0].q[1] = 1;
  jbodies[0].q[2] = 1;
  Cells cells(4);
  Cell * Cj = &cells[0];
  Cj->X[0] = 3;
  Cj->X[1] = 1;
  Cj->X[2] = 1;
  Cj->R = 1;
  Cj->body = &jbodies[0];
  Cj->numBodies = jbodies.size();
  Cj->M.resize(3);
  for(int ind=0; ind<3; ind++) Cj->M[ind].resize(NTERM, 0.0);
  P2M(Cj);

  //! M2M
  Cell * CJ = &cells[1];
  CJ->child = Cj;
  CJ->numChilds = 1;
  CJ->X[0] = 4;
  CJ->X[1] = 0;
  CJ->X[2] = 0;
  CJ->R = 2;
  CJ->M.resize(3);
  for(int ind=0; ind<3; ind++) CJ->M[ind].resize(NTERM, 0.0);
  M2M(CJ);

  //! M2L
  Cell * CI = &cells[2];
  CI->X[0] = -4;
  CI->X[1] = 0;
  CI->X[2] = 0;
  CI->R = 2;
  CI->L.resize(3);
  for(int ind=0; ind<3; ind++) CI->L[ind].resize(NTERM, 0.0);
  M2L(CI, CJ);

  //! L2L
  Cell * Ci = &cells[3];
  CI->child = Ci;
  CI->numChilds = 1;
  Ci->X[0] = -3;
  Ci->X[1] = 1;
  Ci->X[2] = 1;
  Ci->R = 1;
  Ci->L.resize(3);
  for(int ind=0; ind<3; ind++) Ci->L[ind].resize(NTERM, 0.0);
  L2L(CI);

  //! L2P
  Bodies bodies(1);
  bodies[0].X[0] = -2;
  bodies[0].X[1] = 2;
  bodies[0].X[2] = 2;
  bodies[0].q[0] = 1;
  bodies[0].q[1] = 1;
  bodies[0].q[2] = 1;
  bodies[0].p = 0;
  bodies[0].J = 0;
  Ci->body = &bodies[0];
  Ci->numBodies = bodies.size();
  L2P(Ci);

  //! P2P
  Bodies bodies2(1);
  for (size_t b=0; b<bodies2.size(); b++) {
    bodies2[b] = bodies[b];
    bodies2[b].p = 0;
    bodies2[b].J = 0;
  }
  Cj->numBodies = jbodies.size();
  Ci->numBodies = bodies2.size();
  Ci->body = &bodies2[0];
  P2P(Ci, Cj);

  //! Verify results
  Verify verify;
  vec3 pDif = verify.getDifScalar(bodies, bodies2);
  vec3 pNrm = verify.getNrmScalar(bodies2);
  vec3 pRel = 0;
  for(int ind=0; ind<3; ind++) pRel[ind] = std::sqrt(pDif[ind]/pNrm[ind]);
  vec3 FDif = verify.getDifVector(bodies, bodies2);
  vec3 FNrm = verify.getNrmVector(bodies2);
  vec3 FRel = 0;
  for(int ind=0; ind<3; ind++) FRel[ind] = std::sqrt(FDif[ind]/FNrm[ind]);
  std::cout << "Rel. L2 Error (p) : " << pRel << "\n";
  std::cout << "Rel. L2 Error (F) : " << FRel << "\n";
  return 0;
}
