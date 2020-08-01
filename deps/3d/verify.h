#ifndef verify_h
#define verify_h
#include <fstream>
#include "exafmm.h"

namespace exafmm {
  class Verify {
  private:
    const char * path;

  public:
    Verify(const char * _path="./") : path(_path) {}

    vec3 getSumScalar(const Bodies & bodies) {
      vec3 v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++) v[ind] += bodies[b].p[ind] * bodies[b].q[ind];
      }
      return v;
    }

    vec3 getNrmScalar(const Bodies & bodies) {
      vec3 v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++) v[ind] += bodies[b].p[ind] * bodies[b].p[ind];
      }
      return v;
    }

    vec3 getDifScalar(const Bodies & bodies, const Bodies & bodies2) {
      vec3 v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++){
          v[ind] += (bodies[b].p[ind] - bodies2[b].p[ind]) * (bodies[b].p[ind] - bodies2[b].p[ind]);
        }
      }
      return v;
    }

    vec3 getRelScalar(const Bodies & bodies, const Bodies & bodies2) {
      vec3 v = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++){
          v += (bodies[b].p[ind] - bodies2[b].p[ind]) * (bodies[b].p[ind] - bodies2[b].p[ind])
          / (bodies2[b].p[ind] * bodies2[b].p[ind]);
        }
      }
      return v;
    }

    vec3 getNrmVector(const Bodies & bodies) {
      vec3 v = 0;
      vec3 aux1 = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++){
          for(int ind2=0; ind2<3; ind2++) aux1[ind2] = bodies[b].J[ind*3 + ind2];
          v[ind] += norm(aux1);
        }
      }
      return v;
    }

    vec3 getDifVector(const Bodies & bodies, const Bodies & bodies2) {
      vec3 v = 0;
      vec3 aux1 = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++){
          for(int ind2=0; ind2<3; ind2++){
            aux1[ind2] = bodies[b].J[ind*3 + ind2] - bodies2[b].J[ind*3 + ind2];
          }
          v[ind] += norm(aux1);
        }
      }
      return v;
    }

    vec3 getRelVector(const Bodies & bodies, const Bodies & bodies2) {
      vec3 v = 0;
      vec3 aux1 = 0;
      vec3 aux2 = 0;
      for (size_t b=0; b<bodies.size(); b++) {
        for(int ind=0; ind<3; ind++){
          for(int ind2=0; ind2<3; ind2++){
            aux1[ind2] = bodies[b].J[ind*3 + ind2] - bodies2[b].J[ind*3 + ind2];
            aux2[ind2] = bodies2[b].J[ind*3 + ind2];
          }
          v += norm(aux1) / norm(aux2);
        }
      }
      return v;
    }
  };
}
#endif
