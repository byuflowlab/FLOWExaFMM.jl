# FLOWExaFMM
Julia wrapper of ExaFMM with modifications for a vortex solver. Developed and tested in Julia v1.6.1.

First commit in this repo is a clone of https://github.com/EdoAlvarezR/MyExaFMM at 20c6603127dea2d3f2d8200917653a708e0e29f5, which was initiated as a clone of the FLOWVPM_integration branch from a [previous fork of ExaFMM](https://github.com/EdoAlvarezR/exafmm/tree/FLOWVPM_integration) on Aug 10th, 2018.


# Compilation Instructions
To compile the C++ part of this code, first modify the user-provided paths in `build.sh` and then run the build script as `sh build.sh`.

* If running into issue, test that you can correctly compile C++ code wrapped for Julia following these instructions: [LINK](https://nbviewer.jupyter.org/url/edoalvar2.groups.et.byu.net/LabNotebook/202008/FLOWVPMSetupFinal.ipynb)
* Mac user may also need to take a look at these instructions: [LINK](https://github.com/byuflowlab/FLOWUnsteady/issues/26)
* Instructions for Fulton supercomputer: [LINK](https://nbviewer.jupyter.org/url/edoalvar2.groups.et.byu.net/LabNotebook/202108/FLOWVPMSuperComputer.ipynb)

# Automated Compilation (for use with BYU's supercomputer)

Here is an example of an input script I use to send jobs to the supercomputer using `sbatch`.

```
#!/bin/bash

#SBATCH --time=02:00:00   # walltime
#SBATCH --ntasks=24   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --mem-per-cpu=1024M   # memory per CPU core
#SBATCH -J "/fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_01_06.sh"   # job name


# Set the max number of threads to use for programs using OpenMP. Should be <= ppn. Does nothing if the program doesn't use OpenMP.
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load julia/1.6

# build fmm
source /fslhome/rander39/julia/dev/FLOWExaFMM/build_tmp.sh


# include julia files
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_01.jl
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_02.jl
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_03.jl
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_04.jl
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_05.jl
julia /fslhome/rander39/julia/dev/FLOWGround/scripts/20220401/220401_06.jl
```

The idea is `build_tmp.sh` sets environment variables which allows `ExaFMM` to point to the binary built on the node allocated for the job. I have not tested this for use across multiple nodes.

# Authorship
ExaFMM (https://github.com/exafmm) was created and developed by Rio Yokota and Lorena Barba and licensed as BSD 3-Clause.
FLOWExaFMM is a version of ExaFMM with modifications by Eduardo J. Alvarez (Edo.AlvarezR@gmail.com) at BYU's FLOW Lab, to be used in a viscous vortex particle code, and wrapped for the Julia language.

FLOWExaFMM preserves the same BSD 3-Clause license of the original code.
