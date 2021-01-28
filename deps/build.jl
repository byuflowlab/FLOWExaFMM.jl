##############################################
# build script for FLOWExaFMM
# Author: Ryan Anderson
# Date: 9 Sept 2020
##############################################

import CxxWrap

#=
NOTE: for Mac users, the following environment variables must be set before `FLOWExaFMM` will build properly. This script will attempt to find and set them, but some manual configuration may be required.

```bash
export CMAKE_C_COMPILER_JLENV="/usr/local/Cellar/llvm/10.0.1/bin/clang"
export CMAKE_CXX_COMPILER_JLENV="/usr/local/Cellar/llvm/10.0.1/bin/clang++"
export OPENMP_LIBRARIES_JLENV="/usr/local/Cellar/llvm/10.0.1/lib"
export OPENMP_INCLUDES_JLENV="/usr/local/Cellar/llvm/10.0.1/include"
export JlCxx_DIR=/Users/randerson/.julia/artifacts/6017255205dc4fbf4d962903a855a0c631f092dc/
```
=#

# find package directory
println("joinpath(@__DIR__, '..') = ",joinpath(@__DIR__, ".."))
packagedir = normpath(joinpath(@__DIR__, ".."))

if Sys.isapple()
	# check for environment variables
	if isdir("/usr/local/Cellar/llvm/")
		llvm_dir = "/usr/local/Cellar/llvm/"
	else
		throw("Homebrew's LLVM installation not found; run `brew install llvm` and then build `FLOWExaFMM`")
	end
	useversion = ""
    for thisversion in readdir(llvm_dir)
        println("thisversion = ", thisversion)
		global useversion = joinpath(llvm_dir,thisversion)
    end
    println("useversion = ",useversion)
	CMAKE_C_COMPILER_JLENV = joinpath(useversion,"bin/clang")
	CMAKE_CXX_COMPILER_JLENV = joinpath(useversion,"bin/clang++")
	OPENMP_LIBRARIES_JLENV = joinpath(useversion,"lib")
	OPENMP_INCLUDES_JLENV = joinpath(useversion,"include")
    if !(isfile(CMAKE_C_COMPILER_JLENV) && isfile(CMAKE_CXX_COMPILER_JLENV) && isdir(OPENMP_LIBRARIES_JLENV) && isdir(OPENMP_INCLUDES_JLENV))
        throw("The following environment variables could not be set:\n\tCMAKE_C_COMPILER_JLENV\n\tCMAKE_CXX_COMPILER_JLENV\n\tOPENMP_LIBRARIES_JLENV\n\tOPENMP_INCLUDES_JLENV\n\nSet these manually and rebuild `FLOWExaFMM`")
    end
	# set JLCxx prefix path variable
	JlCxx_DIR = CxxWrap.prefix_path()

	# write to file
	buildcontents = open(joinpath(packagedir, "build.sh"), "r") do io
		return read(io, String)
    end
    build_cutoff_i = findfirst(isequal('@'), buildcontents)
    buildcontents = buildcontents[build_cutoff_i+1:end]
    open(joinpath(packagedir, "build.sh"), "w") do io
        write(io, "# set environment variables to use LLVM compilers and libraries\n")
        write(io, "export CMAKE_C_COMPILER_JLENV=$(CMAKE_C_COMPILER_JLENV)\n")
        write(io, "export CMAKE_CXX_COMPILER_JLENV=$(CMAKE_CXX_COMPILER_JLENV)\n")
        write(io, "export OPENMP_LIBRARIES_JLENV=$(OPENMP_LIBRARIES_JLENV)\n")
        write(io, "export OPENMP_INCLUDES_JLENV=$(OPENMP_INCLUDES_JLENV)\n\n")
        write(io, "# set environment variables for use with CxxWrap\n")
        write(io, "export JlCxx_DIR=$(JlCxx_DIR)\n\n# @")
        write(io, buildcontents)
    end

    makecontents = open(joinpath(packagedir, "make.sh"), "r") do io
		return read(io, String)
    end
    make_cutoff_i = findfirst(isequal('@'), makecontents)
    makecontents = makecontents[make_cutoff_i+1:end]
    open(joinpath(packagedir, "make.sh"), "w") do io
        write(io, "# set environment variables to use LLVM compilers and libraries\n")
        write(io, "export CMAKE_C_COMPILER_JLENV=$(CMAKE_C_COMPILER_JLENV)\n")
        write(io, "export CMAKE_CXX_COMPILER_JLENV=$(CMAKE_CXX_COMPILER_JLENV)\n")
        write(io, "export OPENMP_LIBRARIES_JLENV=$(OPENMP_LIBRARIES_JLENV)\n")
        write(io, "export OPENMP_INCLUDES_JLENV=$(OPENMP_INCLUDES_JLENV)\n\n")
        write(io, "# set environment variables for use with CxxWrap\n")
        write(io, "export JlCxx_DIR=$(JlCxx_DIR)\n\n# @")
        write(io, buildcontents)
    end

    # println("CMAKE_C_COMPILER_JLENV = ", CMAKE_C_COMPILER_JLENV)
    # println("CMAKE_CXX_COMPILER_JLENV = ", CMAKE_CXX_COMPILER_JLENV)
    # println("OPENMP_LIBRARIES_JLENV = ", OPENMP_LIBRARIES_JLENV)
    # println("OPENMP_INCLUDES_JLENV = ", OPENMP_INCLUDES_JLENV)
	# if isfile(CMAKE_C_COMPILER_JLENV) && isfile(CMAKE_CXX_COMPILER_JLENV) && isdir(OPENMP_LIBRARIES_JLENV) && isdir(OPENMP_INCLUDES_JLENV)
	#     buildcontents = open(joinpath(packagedir, "build.sh"), "r") do io
	# 	    return read(io, String)
	#     end
	#     open(joinpath(packagedir, "build.sh"), "w") do io
	# 	    write(io, "# set environment variables to use LLVM compilers and libraries\n")
	# 	    write(io, "export CMAKE_C_COMPILER_JLENV=$(CMAKE_C_COMPILER_JLENV)\n")
	# 	    write(io, "export CMAKE_CXX_COMPILER_JLENV=$(CMAKE_CXX_COMPILER_JLENV)\n")
	# 	    write(io, "export OPENMP_LIBRARIES_JLENV=$(OPENMP_LIBRARIES_JLENV)\n")
	# 	    write(io, "export OPENMP_INCLUDES_JLENV=$(OPENMP_INCLUDES_JLENV)\n")
	# 	    write(io, buildcontents)
	#     end
	# else
	# 	println("The following environment variables could not be set:\n\tCMAKE_C_COMPILER_JLENV\n\tCMAKE_CXX_COMPILER_JLENV\n\tOPENMP_LIBRARIES_JLENV\n\tOPENMP_INCLUDES_JLENV\n\nSet these manually and rebuild `FLOWExaFMM`")
	# end

end

buildcommand = joinpath(packagedir, "build.sh")
makecommand = joinpath(packagedir, "make.sh")
# run(`sh -c $buildcommand`)
# run(`sh -c $makecommand`)
println("Finished building FLOWExaFMM")
