#ifndef buildtree_h
#define buildtree_h
#include "exafmm.h"

namespace exafmm {

  //! Calculate cell's characteristic smoothing radius
  void calcSigmas(Cell * Cj){

    Cj->sigma = 0.0;

    // Recursive call
    for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
#pragma omp task untied if(Ci->numBodies > 100)
      calcSigmas(Ci);
    }
#pragma omp taskwait

    // Base case
    if (Cj->numChilds==0){
      float_t weights =  0.0;
      for (Body * B=Cj->body; B!=Cj->body+Cj->numBodies; B++) {
        Cj->sigma += norm(B->q)*B->sigma[0];
        weights += norm(B->q);
      }
      Cj->sigma = Cj->sigma/weights;
    }

    // Calculates this cell's sigma as the average of its populated childs' sigmas
    else{
      int nchilds = 0;
      for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
        if (Ci->numBodies!=0){
          nchilds += 1;
          Cj->sigma += Ci->sigma;
        }
      }
      if (nchilds!=0) Cj->sigma = Cj->sigma/nchilds;
    }
  }

  void calcSigmas(Cells & cells){
    calcSigmas(&cells[0]);
  }

  //! Get bounding box of bodies
  void getBounds(Bodies & bodies, int nb, real_t & R0, vec3 & X0) {
    vec3 Xmin = bodies[0].X;
    vec3 Xmax = bodies[0].X;
    for (int i=0; i<nb; i++) {
      Xmin = min(bodies[i].X, Xmin);
      Xmax = max(bodies[i].X, Xmax);
    }
    X0 = (Xmax + Xmin) / 2;
    R0 = fmax(max(X0-Xmin), max(Xmax-X0));
    R0 *= 1.00001;
  }

  //! Build cells of tree adaptively using a top-down approach based on recursion
  void buildCells(Body * bodies, Body * buffer, int begin, int end, Cell * cell, Cells & cells,
                  const vec3 & X, real_t R, int level=0, bool direction=false) {

    //! Create a tree cell
    cell->body = bodies + begin;
    if(direction) cell->body = buffer + begin;
    cell->numBodies = end - begin;
    cell->numChilds = 0;
    cell->X = X;
    cell->R = R / (1 << level);
    //! If cell is a leaf
    if (end - begin <= NCRIT) {
      if (direction) {
        for (int i=begin; i<end; i++) {
          buffer[i].X = bodies[i].X;
          buffer[i].q = bodies[i].q;
          // buffer[i].vort = bodies[i].vort;
          buffer[i].sigma = bodies[i].sigma;
          buffer[i].vol = bodies[i].vol;
          buffer[i].index = bodies[i].index;
          // buffer[i].p = bodies[i].p;
          if(SGS){
            buffer[i].J = bodies[i].J;
            buffer[i].dJdx1 = bodies[i].dJdx1;
            buffer[i].dJdx2 = bodies[i].dJdx2;
            buffer[i].dJdx3 = bodies[i].dJdx3;
          }
          // buffer[i].pse = bodies[i].pse;
        }
      }
      return;
    }
    //! Count number of bodies in each octant
    int size[8] = {0,0,0,0,0,0,0,0};
    vec3 x;
    for (int i=begin; i<end; i++) {
      x = bodies[i].X;
      int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
      size[octant]++;
    }
    //! Exclusive scan to get offsets
    int offset = begin;
    int offsets[8], counter[8];
    for (int i=0; i<8; i++) {
      offsets[i] = offset;
      offset += size[i];
      if (size[i]) cell->numChilds++;
    }
    //! Sort bodies by octant
    for (int i=0; i<8; i++) counter[i] = offsets[i];
    for (int i=begin; i<end; i++) {
      x = bodies[i].X;
      int octant = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
      buffer[counter[octant]].X = bodies[i].X;
      buffer[counter[octant]].q = bodies[i].q;
      // buffer[counter[octant]].vort = bodies[i].vort;
      buffer[counter[octant]].sigma = bodies[i].sigma;
      buffer[counter[octant]].index = bodies[i].index;
      buffer[counter[octant]].vol = bodies[i].vol;
      // buffer[counter[octant]].p = bodies[i].p;
      if(SGS){
        buffer[counter[octant]].J = bodies[i].J;
        buffer[counter[octant]].dJdx1 = bodies[i].dJdx1;
        buffer[counter[octant]].dJdx2 = bodies[i].dJdx2;
        buffer[counter[octant]].dJdx3 = bodies[i].dJdx3;
      }
      // buffer[counter[octant]].pse = bodies[i].pse;
      counter[octant]++;
    }
    //! Loop over children and recurse
    vec3 Xchild;
    cells.resize(cells.size()+cell->numChilds);
    Cell * child = &cells.back() - cell->numChilds + 1;
    cell->child = child;
    int c = 0;
    for (int i=0; i<8; i++) {
      Xchild = X;
      real_t r = R / (1 << (level + 1));
      for (int d=0; d<3; d++) {
        Xchild[d] += r * (((i & 1 << d) >> d) * 2 - 1);
      }
      if (size[i]) {
        buildCells(buffer, bodies, offsets[i], offsets[i] + size[i],
                   &child[c], cells, Xchild, R, level+1, !direction);
        c++;
      }
    }
  }

  Cells buildTree(Bodies & bodies, int nb) {
    real_t R0;
    vec3 X0;
    getBounds(bodies, nb, R0, X0);
    Bodies buffer = bodies;
    Cells cells(1);
    cells.reserve(nb);

    buildCells(&bodies[0], &buffer[0], 0, nb, &cells[0], cells, X0, R0);

    calcSigmas(cells);

    return cells;
  }
}
#endif
