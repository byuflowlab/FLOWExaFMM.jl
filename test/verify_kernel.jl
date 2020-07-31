vpm_path = "/home/user/Dropbox/FLOWResearch/MyCodes/MyVPM/"
# vpm_path = "../../MyVPM20180814/"
include(joinpath(vpm_path*"src/MyVPM.jl"))
vpm = MyVPM

module_path = splitdir(@__FILE__)[1]

using PyPlot


function test_kernels(; n_particles::Int64=200, probe_i::Int64=1,
                        # EXAFMM OPTIONS
                        P=vpm.exafmm_P, NCRIT=vpm.exafmm_NCRIT,
                        THETA=vpm.exafmm_THETA, PHI=vpm.exafmm_PHI,
                        # TEST SETUP
                        box_length=1.0, Gamma_low=zeros(3), Gamma_up=ones(3),
                        sigma_low=0.02, sigma_up=0.2,
                        # TEST OPTIONS
                        verbose=true, prompt=false, disp_plot=true,
                        kernels_to_test=["Direct Gaussian", "FMM Gaussian"],
                        # kernels_to_test=["Direct Singular", "FMM Singular"],
                        # kernels_to_test=["Direct Winckelmann's", "FMM Winckelmann's",
                        #                     "FMM Winckelmann's noH"],
                        refkernel=1)

  vpm.set_P(Int32(P))
  vpm.set_NCRIT(Int32(NCRIT))
  vpm.set_THETA(Float64(THETA))
  vpm.set_PHI(Float64(PHI))

  # PARAMETERS
  N = n_particles             # Number of particles
  X_low = -box_length*ones(3) # Lower bound of particle position
  X_up = -X_low               # Upper bound of particle position

  # Gamma_low = zeros(3)      # Lower bound of vectorial circulation
  # Gamma_up = ones(3)        # Upper bound of vectorial circulation

  # sigma_low = 0.02            # Lower bound of smoothing radius
  # sigma_up = 0.2              # Upper bound of smoothing radius

  # Generates the particle field
  Vinf(X, t) = zeros(3)
  pfield = vpm.ParticleField(N, Vinf, nothing, "ExaFMM")
  for i in 1:N
    particle = vcat(
                      X_low+rand(3).*(X_up-X_low),
                      Gamma_low+rand(3).*(Gamma_up-Gamma_low),
                      sigma_low+rand()*(sigma_up-sigma_low),
                      0
                    )
    vpm.addparticle(pfield, particle)
  end

  # Probe
  probe_X = vpm.get_x(pfield, probe_i)
  probe_Gamma = vpm.get_gamma(pfield, probe_i)

  # Reference kernel
  avail_kernels = Dict(
        "Direct Singular" =>  Dict(
                                    "label"     => "Direct Singular",
                                    "exafmm"    => false,
                                    "options"   => ["direct"]
                                  ),
        "Direct Winckelmann's" =>  Dict(
                                    "label"     => "Direct Winckelmann's",
                                    "exafmm"    => false,
                                    "options"   => ["directblobwinckelmann"]
                                  ),
        "Direct Gaussian" =>  Dict(
                                    "label"     => "Direct Gaussian",
                                    "exafmm"    => false,
                                    "options"   => ["directblob"]
                                  ),
        "FMM Singular"    =>  Dict(
                                    "label"     => "FMM Singular",
                                    "exafmm"    => true,
                                    "options"   => ["ExaFMM", Int32(1), Int32(1) ]
                                  ),
        "FMM Winckelmann's"    =>  Dict(
                                    "label"     => "FMM Winckelmann's",
                                    "exafmm"    => true,
                                    "options"   => ["ExaFMM", Int32(3), Int32(1) ]
                                  ),
        "FMM Winckelmann's noH"    =>  Dict(
                                    "label"     => "FMM Winckelmann's noH",
                                    "exafmm"    => true,
                                    "options"   => ["ExaFMM", Int32(3), Int32(0) ]
                                  ),
        "FMM Gaussian"    =>  Dict(
                                    "label"     => "FMM Gaussian",
                                    "exafmm"    => true,
                                    "options"   => ["ExaFMM", Int32(5), Int32(1) ]
                                  ),
        )

  kernels = []

  # Kernels to test
  for k_name in kernels_to_test
    if !(k_name in keys(avail_kernels))
      error("Invalid kernel $(k_name). Available kernels are: $(keys(avail_kernels))")
    else
      push!(kernels, avail_kernels[k_name])
    end
  end

  # Reorganizes them to place reference kernel first
  kernels = vcat(kernels[refkernel], kernels[1:refkernel-1],
                                                      kernels[refkernel+1:end])


  if disp_plot
    fig = figure("verif_error", figsize=(7*2, 5*1))
    suptitle("Field error relative to $(kernels[1]["label"])")
  end


  if verbose
    println("==============================================================")
    println("\tKERNEL COMPARED TO $(kernels[1]["label"])")
    println("==============================================================")
    println("Number of particles: $N")
    println("Probe index: $probe_i")
    println("Probe position: $probe_X")
    println("Multipole order: $(vpm.exafmm_P)")
    println("Multipole NCRIT: $(vpm.exafmm_NCRIT)")
    println("Multipole THETA: $(vpm.exafmm_THETA)")
    println("Multipole PHI: $(vpm.exafmm_PHI)")
    println("==============================================================")
  end

  Uref = nothing
  stchref = nothing
  for (i,kernel) in enumerate(kernels)
    if verbose
      println("**************************************************************")
      println("************ "*kernel["label"]*"*"^(63-14-length(kernel["label"])))
      println("**************************************************************")
    end

    # Kernel setup
    if kernel["exafmm"]
      vpm.set_P2PTYPE(kernel["options"][2])
      vpm.set_L2PTYPE(kernel["options"][3])
    end

    if verbose; println("\tWarm-up run..."); end;
    @time vpm.conv(pfield, kernel["options"][1])

    if verbose; println("\tTiming run..."); end;
    @time M, _ = vpm.conv(pfield, kernel["options"][1])

    U = M[:, 4:6]
    stch = M[:, 1:3]

    # Case of reference kernel
    if i==1
      Uref = U
      stchref = stch

      if verbose
        println("\t********************************************")
        println("\t\tReference U:\t$(Uref[probe_i,:])")
        println("\t\tReference stretch:\t$(stchref[probe_i,:])")
        println("\t********************************************")
      end
    end


    # Calculates error
    Uerr = U - Uref
    stcherr = stch - stchref
    Urms = sqrt(mean([norm(Uerr[i,:])^2 for i in 1:N]))
    stchrms = sqrt(mean([norm(stcherr[i,:])^2 for i in 1:N]))


    if disp_plot && i!=1

      subplot(121)
      # y = [norm(Uerr[i,:])/norm(Uref[i,:])*100 for i in 1:N]
      y = [ mean(abs.(Uerr[i,:]./Uref[i,:])) for i in 1:N]
      plot(collect(1:N), y, "o", alpha=0.5, label=kernel["label"])

      subplot(122)
      # y = [norm(stcherr[i,:])/norm(stchref[i,:])*100 for i in 1:N]
      y = [ mean(abs.(stcherr[i,:]./stchref[i,:])) for i in 1:N]
      plot(collect(1:N), y, "o", alpha=0.5, label=kernel["label"])
    end

    if verbose
      println("\t\tProbe U:\t$(U[probe_i,:])")
      println("\t\tProbe stretch:\t$(stch[probe_i,:])")

      Uerr = norm(U[probe_i,:] - Uref[probe_i,:])/norm(Uref[probe_i,:])
      stcherr = norm(stch[probe_i,:] - stchref[probe_i,:])/norm(stchref[probe_i,:])
      println("\t\tU probe error:\t\t$(@sprintf("%0.2e", Uerr*100 ))%")
      println("\t\tStretch probe error:\t$(@sprintf("%0.2e", stcherr*100 ))%")

      println("\t\tU RMS:\t\t\t$(@sprintf("%0.2e", Urms ))")
      println("\t\tStretch RMS:\t\t$(@sprintf("%0.2e", stchrms ))")
    end

  end

  if disp_plot
    subplot(121)
    xlabel("Particle index")
    ylabel(L"$\Vert {\bf u} \Vert$ error %")
    grid(true, color="0.8", linestyle="--")
    ax = gca()
    ax[:set_yscale]("log")

    subplot(122)
    xlabel("Particle index")
    ylabel(L"$\Vert \left( \Gamma \cdot \nabla \right){\bf u} \Vert$ error %")
    grid(true, color="0.8", linestyle="--")
    legend(loc="upper left", bbox_to_anchor=(1, 0.75))
    ax = gca()
    ax[:set_yscale]("log")
  end
end
