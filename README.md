# FLOWExaFMM
Julia wrapper of ExaFMM with modifications for a vortex solver. Developed and tested in Julia v1.6.

[ExaFMM](https://joss.theoj.org/papers/10.21105/joss.03145) is a fast multipole method for exascale computing developed by Lorena Barba and Rio Yokota

> First commit in this repo is a clone of https://github.com/EdoAlvarezR/MyExaFMM at 20c6603127dea2d3f2d8200917653a708e0e29f5, which was initiated as a clone of the FLOWVPM_integration branch from a [previous fork of ExaFMM](https://github.com/EdoAlvarezR/exafmm/tree/FLOWVPM_integration) on Aug 10th, 2018.


# Compilation Instructions

To compile the C++ part of this package, run the build script as `sh build.sh` (or `sh build_macos.sh` for MacOs), which should create a binary file under `src/`. Try the following if encountering errors:

* Test that you can correctly compile C++ code wrapped for Julia following these instructions: [LINK](https://nbviewer.org/github/byuflowlab/FLOWVPM.jl/blob/master/docs/installation-linux.ipynb)
* Mac user may also need to take a look at these instructions: [LINK](https://github.com/byuflowlab/FLOWUnsteady/issues/26)
* Instructions for BYU's Fulton supercomputer: [LINK](https://nbviewer.jupyter.org/url/edoalvar2.groups.et.byu.net/LabNotebook/202108/FLOWVPMSuperComputer.ipynb)


# Authorship
ExaFMM (https://github.com/exafmm) was created and developed by Rio Yokota and Lorena Barba and licensed as BSD 3-Clause.
FLOWExaFMM is a version of ExaFMM with modifications by Eduardo J. Alvarez (Edo.AlvarezR@gmail.com) at BYU's FLOW Lab, to be used in a viscous vortex particle code, and wrapped for the Julia language.

FLOWExaFMM preserves the same BSD 3-Clause license of the original code.
