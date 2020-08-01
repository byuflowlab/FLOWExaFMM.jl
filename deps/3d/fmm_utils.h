
// Generate bodies out of array data
Bodies genBodies(double data[][_n_param], int np){
  Bodies bodies(np);
  for(int i=0; i<np; ++i){
    bodies[i].X[0] = data[i][_x_i-1];
    bodies[i].X[1] = data[i][_y_i-1];
    bodies[i].X[2] = data[i][_z_i-1];
    bodies[i].q[0] = data[i][_gamma_x_i-1];
    bodies[i].q[1] = data[i][_gamma_y_i-1];
    bodies[i].q[2] = data[i][_gamma_z_i-1];
    bodies[i].sigma= data[i][_blob_sigma_i-1];
    bodies[i].vol= data[i][_vol_i-1];
    bodies[i].index= i+1;
  }
  return bodies;
}

Bodies genBodies(double **data, int np){
  Bodies bodies(np);
  for(int i=0; i<np; ++i){
    bodies[i].X[0] = data[i][_x_i-1];
    bodies[i].X[1] = data[i][_y_i-1];
    bodies[i].X[2] = data[i][_z_i-1];
    bodies[i].q[0] = data[i][_gamma_x_i-1];
    bodies[i].q[1] = data[i][_gamma_y_i-1];
    bodies[i].q[2] = data[i][_gamma_z_i-1];
    bodies[i].sigma= data[i][_blob_sigma_i-1];
    bodies[i].vol= data[i][_vol_i-1];
    bodies[i].index= i+1;
  }
  return bodies;
}


// Reads a .pf file
Bodies readBodies(const char file_name[], bool verbose=false){
  std::string line;
  std::ifstream f(read_file, std::ios::in);

  // ERROR CASE
  if (f.is_open()==false){
    std::cout << "Unable to open file";
    exit(1);
  }

  // ---- READS THE HEADER
  for(int i=1; i<=header_length; ++i){
    getline(f,line);
    // std::cout << line << '\n';
    if(i==3) np = atoi(line.c_str());
    else if(i==6) nparam = atoi(line.c_str());
  }
  if(verbose){
    std::cout << "Number of particles: " << np << '\n';
    std::cout << "Number of parameters: " << nparam << '\n';
  }

  if( _n_param!=nparam ){
    std::cout << "********* WARNING: _n_param!=nparam" << '\n';
  }

  // ---- READS THE FILE
  double data[np][_n_param];     // Parameters of each particle
  size_t pos;                  // Parcing of each line
  std::string elem;            // Ditto
  int j;                       // Ditto
  // Iterates over each line
  for(int i=0; i<np; ++i){
    // std::cout << "\nLINE "<< i << '\n' << '\t';
    getline(f,line);
    line += delimiter;

    pos = 0;
    j = 0;
    // Iterates over each element in the line
    while ((pos = line.find(delimiter)) != std::string::npos){
        elem = line.substr(0, pos);
        data[i][j++] = strtod(elem.c_str(), NULL);
        // std::cout << elem << '\t';
        line.erase(0, pos + 1);
    }
    // Error case
    if(j!=(nparam)){
      std::cout << "Logic error 1!";
      exit(1);
    }
  }
  f.close();

  return genBodies(data, np);
}


void closeConnectivity_ilist(Cell * Ctarg, Cell * Csour, std::vector<Cell*> & ilist){
  real_t r2 = norm(Ctarg->X-Csour->X);
  real_t aux1 = (sqrt(r2) - (Ctarg->R+Csour->R));
  // Case that the clusters are sufficiently close for overlap
  if ( aux1<=0 || (Ctarg->sigma+Csour->sigma)/2 / aux1  >= LAMBDA || (Ctarg->R+Csour->R)*(Ctarg->R+Csour->R)/r2 > THETA*THETA ){
// if (
//       ( (r2 - (Ci.R+Cj.R)**2) > 0 && ((Ci.sigma+Cj.sigma)/2/LAMBDA + (Ci.R+Cj.R))**2 < r2)
//       ||
//       ( (r2 - (Ci.R+Cj.R)**2) <= 0 && ((Ci.sigma+Cj.sigma)/2/LAMBDA + (Ci.R+Cj.R))**2 >= r2)
//     ){

    if (Csour->numChilds == 0){
      ilist.push_back(Csour);
    }
    else{
      for (Cell * Ci=Csour->child; Ci!=Csour->child+Csour->numChilds; Ci++){
        closeConnectivity_ilist(Ctarg, Ci, ilist);
      }
    }
  }
}


/*
  Creates a list of bodies closing connectivity issues where particles are
separating and the coreoverlap is not enough. Add this bodies to the particle
field to take care of Lagrangian distortion.
*/
void closeConnectivity(Cell * Croot, Cell * Cj, Bodies & bodies){
  // TODO: Construct Bodies as a list push_back() instead

  bool flag;
  bool flagbd;
  // vec3 qdir;
  real_t R2;
  real_t h;
  Body Btarg;

  // Recursive call
  for (Cell * Ci=Cj->child; Ci!=Cj->child+Cj->numChilds; Ci++) {
// #pragma omp task untied if(Ci->numBodies > 100)
    closeConnectivity(Croot, Ci, bodies);
  }
// #pragma omp taskwait

  // Base case
  if (Cj->numChilds==0){
    Body * Bj = Cj->body;

    std::vector<Cell*> ilist(0);                 // interaction list of this cluster
    closeConnectivity_ilist(Cj, Croot, ilist);

    // Iterates over every body in the cluster checking for connectivity
    if(Cj->numBodies>0) Btarg = Bj[0];
    for (int j=0; j<Cj->numBodies; j++) {

      flag = false;       // Indicates whether connectivity is satisfied
      flagbd = true;      // Indicates whether the body is at the boundary of computational domain
      vec3 qdir = Btarg.q *1/sqrt(norm(Btarg.q));

      // Checks against every body in the same cluster
      // Body * Bi = Cj->body;
      // for (int i=0; !flag && i<Cj->numBodies; i++) {
      //   R2 = norm(Bi[i].X - Btarg.X);
      //   h = dot(qdir, Bi[i].X - Btarg.X);
      //   if (R2 != 0 && h >= 0){                                                 // Ahead condition
      //     flag = (Btarg.sigma+Bi[i].sigma)/2*(Btarg.sigma+Bi[i].sigma)/2 / R2 > LAMBDA*LAMBDA;                    // Core overlap condition
      //     flagbd = flagbd && norm((Bi[i].X - Btarg.X)-qdir*h)/h/h < 0.333;              // Not C. dom. condition: 30deg angle
      //   }
      //
      //   if(Btarg.index==150 && Bi[i].index==151) std::cout << "\t\t******************* TRUE RABBIT1!\t"<<flag<<"\n";
      // }

      // Checks bodies in nearby clusters (including itself)
        for (size_t i=0; !flag && i<ilist.size(); i++) {
          Body * Bk = ilist[i]->body;
          for (int k=0; !flag && k<ilist[i]->numBodies; k++) {
            R2 = norm(Bk[k].X - Btarg.X);
            h = dot(qdir, Bk[k].X - Btarg.X);
            if (R2 != 0 && h >= 0){                                                 // Ahead condition
              flag = (Btarg.sigma+Bk[k].sigma)/2*(Btarg.sigma+Bk[k].sigma)/2 / R2 > LAMBDA*LAMBDA;                    // Core overlap condition
              flagbd = flagbd && norm((Bk[k].X - Btarg.X)-qdir*h)/h/h < 0.333;              // Not C. dom. condition: 30deg angle
            }

            // if(Btarg.index==150) std::cout << "\t\ti="<<i<<"\tBk[k].index="<<Bk[k].index<<"\n";
            // if(Btarg.index==150 && Bk[k].index==151) std::cout << "\t\t******************* TRUE RABBIT2!\t"<<flag<<"\n";
          }
        }

      // Checks bodies to be added
      for (size_t i=0; !flag && i<bodies.size(); i++) {
        R2 = norm(bodies[i].X - Btarg.X);
        h = dot(qdir, bodies[i].X - Btarg.X);
        if (R2 != 0 && h >= 0){                                                 // Ahead condition
          flag = (Btarg.sigma+bodies[i].sigma)/2*(Btarg.sigma+bodies[i].sigma)/2 / R2 > LAMBDA*LAMBDA;                    // Core overlap condition
          flagbd = flagbd && norm((bodies[i].X - Btarg.X)-qdir*h)/h/h < 0.333;              // Not C. dom. condition: 30deg angle
        }
      }

    // if (!flag) std::cout << "\t\t******************* TRUE RABBIT!\n";

      // It adds a body if connectivity wasn't satisfied and is not in comp. dom. boundary
      if (!flag && !flagbd){
        Body body;
        body.X = Btarg.X + qdir*Btarg.sigma/LAMBDA;
        // Since placing the new particle right at the boundary where it
        // satisfies lambda may lead to too much overlap with a nearby particle,
        // here I'm placing it half way instead
        // body.X = Btarg.X + qdir*Btarg.sigma/LAMBDA/2;
        // body.q[0] = 0;
        // body.q[1] = 0;
        // body.q[2] = 0;
        // Here it gives it a small strength to point in the direction of connectivity
        // body.q = qdir*EPS;
        body.q = qdir;
        body.sigma= Btarg.sigma;
        body.vol= Btarg.vol;
        body.index= -1;
        bodies.push_back(body);

        // Loops again checking whether the new body satisfies connectivity
        Btarg = body;
        j--;
      }
      else if(j+1<Cj->numBodies){
        Btarg = Bj[j+1];
      }
    }
  }
}
Bodies closeConnectivity(Cells & cells){
  Bodies bodies(0);
  int numbodies = 0;
  closeConnectivity(&cells[0], &cells[0], bodies);

  initTarget(bodies);
  return bodies;
}
