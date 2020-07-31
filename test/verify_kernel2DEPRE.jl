VPM_PATH = "/home/user/Dropbox/FLOWResearch/MyCodes/MyVPM/src/"
# VPM_PATH = "MyVPM/"
include(VPM_PATH*"MyVPM.jl")
vpm = MyVPM

using PyPlot
pygui(true)

function run(n_particles; probe_i=1, verbose=true, plot_=false,
                    P::Int32=vpm.exafmm_P, NCRIT::Int32=vpm.exafmm_NCRIT,
                    THETA::Float64=vpm.exafmm_THETA, fmmverbose=false,
                    error="velocity", sigmablob=0.1,
                    n_probe_particles=0,
                    dirblob=true, dirsing=false)
  vpm.set_P(P)
  vpm.set_NCRIT(NCRIT)
  vpm.set_THETA(THETA)
  vpm.set_VERBOSE(fmmverbose)

  # METHODS
  # dirblob = true
  # dirsing = false
  exafmm = false
  exafmmfd = true
  exafmmsingANAcsda = false
  exafmmregANAcsda = true
  exafmmcsda = false

  # PARAMETERS
  len = 50                    # Length of box edge
  C = [100, -50, 25]          # Center of box
  l = C - len*[0.5,0.5,0.5]   # Lower bound for spawning particles
  u = C + len*[0.5,0.5,0.5]   # Upper bound for spawning particles

  qlen = 1                    # Dispersion of charge
  qC = [1, 1, 1]*qlen*0.75    # Center of charge
  ql = qC - qlen*[0.5,0.5,0.5]# Lower bound of charges
  qu = qC + qlen*[0.5,0.5,0.5]# Upper bound of charges

  # n_probe_particles = 100     # Number of probe particles
  probe_fact = 0.5           # Factor outside of the box to probe

  # n_particles = 100           # Number of particles
  # probe_i = 1                  # Index of the particle that will be the probe

  # Generates the particles
  Q = Array{Float64,1}[]        # Position
  q = Array{Float64,1}[]        # Charge
  for i in 1:n_particles
      this_pos = l + rand(3).*(u-l)
      this_chrg = ql + rand(3).*(qu-ql)
      push!(Q, this_pos)
      push!(q, this_chrg)
  end

  # Generates particle field
  Uinf(X,t ) = [0.0,0.0,0.0]
  p_field = vpm.ParticleField(n_particles+n_probe_particles, Uinf, nothing,
                              "directblob")
  for i in 1:n_particles
    vpm.addparticle(p_field, vcat(Q[i], q[i], sigmablob) )
  end


  # Line of probes
  probe_Us, probe_convs = [], []
  part_Us, part_convs = [], []
  method_labels = []        # Use this to keep track of method orders
  aux1 = (u-l)
  pXmin = l - probe_fact*aux1
  pXmax = u + probe_fact*aux1
  aux2 = norm(pXmax-pXmin)
  aux3 = (pXmax-pXmin)/aux2
  # Generates probes
  for i in 1:n_probe_particles
      this_pos = pXmin + i*(aux2/n_probe_particles)*aux3
      this_chrg = zeros(3)
      push!(Q, this_pos)
      push!(q, this_chrg)
  end
  for i in n_particles+1:n_particles + n_probe_particles
    vpm.addparticle(p_field, vcat(Q[i], q[i], sigmablob) )
  end

  # Individual probe
  probe_X = vpm.get_x(p_field, probe_i)

  if verbose
    println("==============================================================")
    println("Number of particles: $n_particles")
    println("Probe index: $probe_i")
    println("Probe position: $probe_X")
    println("Multipole order: $(vpm.exafmm_P)")
    println("Multipole NCRIT: $(vpm.exafmm_NCRIT)")
    println("Multipole THETA: $(vpm.exafmm_THETA)")
    println("==============================================================")
  end

    Uprobes = [[0.0,0,0], [1.0, 1.0, 1.0]] # Just testing the probe features

    if dirblob
      # ------------ MyVPM Direct Blob
      this_label = "DirectBlob"
      push!(method_labels, this_label)
      if verbose;println("\t***************** $this_label *****************");end;

      if verbose; println("\tTest run..."); end;
      vpm.conv(p_field, "directblob"; Uprobes=Uprobes)

      if verbose; println("\tNow for real..."); end;
      prev_t = time()
      Mdir, _ = vpm.conv(p_field, "directblob"; Uprobes=Uprobes)
      end_t = time()
      tdir = end_t-prev_t

      # Induced velocity
      Udir = Mdir[probe_i,4:6]
      convdir = Mdir[probe_i,1:3]

      # Potential
      potdir = zeros(3)
      for i in 1:n_particles
        if i!=probe_i
          dX = vpm.get_x(p_field,probe_i) - vpm.get_x(p_field,i)
          potdir += vpm.get_gamma(p_field,i) / norm(dX)
        end
      end

      # Stores particle values
      Us, convs = [], []
      for i in 1:n_particles
        push!(Us, Mdir[i,4:6])
        push!(convs, Mdir[i,1:3])
      end
      push!(part_Us, Us)
      push!(part_convs, convs)
      # Stores probe values
      Us, convs = [], []
      for i in n_particles+1:n_particles + n_probe_particles
        push!(Us, Mdir[i,4:6])
        push!(convs, Mdir[i,1:3])
      end
      push!(probe_Us, Us)
      push!(probe_convs, convs)

      if verbose
        println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$potdir")
          println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
          println("\t\t\t$Udir")
          println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
          println("\t\t\t$convdir")
        println("\tTime: $(end_t-prev_t) (s)")
      end
    end


  if dirsing
    # ------------ MyVPM Direct Singular
    this_label = "DirectSing"
    push!(method_labels, this_label)
    if verbose;println("\t***************** $this_label *****************");end;

    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "direct")

    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mdir, _ = vpm.conv(p_field, "direct")
    end_t = time()
    tdir = end_t-prev_t

    # Induced velocity
    Udir = Mdir[probe_i,4:6]
    convdir = Mdir[probe_i,1:3]

    # Potential
    potdir = zeros(3)
    for i in 1:n_particles
      if i!=probe_i
        dX = vpm.get_x(p_field,probe_i) - vpm.get_x(p_field,i)
        potdir += vpm.get_gamma(p_field,i) / norm(dX)
      end
    end

    # Stores particle values
    Us, convs = [], []
    for i in 1:n_particles
      push!(Us, Mdir[i,4:6])
      push!(convs, Mdir[i,1:3])
    end
    push!(part_Us, Us)
    push!(part_convs, convs)
    # Stores probe values
    Us, convs = [], []
    for i in n_particles+1:n_particles + n_probe_particles
      push!(Us, Mdir[i,4:6])
      push!(convs, Mdir[i,1:3])
    end
    push!(probe_Us, Us)
    push!(probe_convs, convs)

    if verbose
      println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$potdir")
        println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$Udir")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convdir")
      println("\tTime: $(end_t-prev_t) (s)")
    end
  end

  if exafmm
    # ---------------- ExaFMM
    this_label = "barebone regExaFMM"
    push!(method_labels, this_label)
    if verbose;println("\t***************** $this_label *****************");end;
    vpm.set_P2PTYPE(Int32(2))
    vpm.set_L2PTYPE(Int32(0))

    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "ExaFMM")
    _, _, _, _, pfield_out = vpm.U_ExaFMM(p_field._p_field, vpm.get_np(p_field))
    potexafmm = pfield_out[probe_i][vpm._n_param+1:vpm._n_param+3]
    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm, _ = vpm.conv(p_field, "ExaFMM")
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    # Stores particle values
    Us, convs = [], []
    for i in 1:n_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(part_Us, Us)
    push!(part_convs, convs)
    # Stores probe values
    Us, convs = [], []
    for i in n_particles+1:n_particles + n_probe_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(probe_Us, Us)
    push!(probe_convs, convs)

    if verbose
      println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$potexafmm")
        println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end
  end


  if exafmmfd
    # ---------------- DEPRICATED ExaFMM
    this_label = "FD ExaFMM"
    push!(method_labels, this_label)
    if verbose;println("\t***************** $this_label *****************");end;
    vpm.set_P2PTYPE(Int32(2))
    vpm.set_L2PTYPE(Int32(0))
    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "FDExaFMM"; Uprobes=Uprobes)

    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm, _ = vpm.conv(p_field, "FDExaFMM"; Uprobes=Uprobes)
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    # Stores particle values
    Us, convs = [], []
    for i in 1:n_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(part_Us, Us)
    push!(part_convs, convs)
    # Stores probe values
    Us, convs = [], []
    for i in n_particles+1:n_particles + n_probe_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(probe_Us, Us)
    push!(probe_convs, convs)

    if verbose
      println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end
  end

  if exafmmsingANAcsda
    # ---------------- ExaFMM singANA-CSDA
    this_label = "singANA-CSDA ExaFMM"
    push!(method_labels, this_label)
    if verbose;println("\t***************** $this_label *****************");end;
    vpm.set_P2PTYPE(Int32(1))
    vpm.set_L2PTYPE(Int32(1))

    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "ExaFMM")
    _, _, _, _, pfield_out = vpm.U_ExaFMM(p_field._p_field, vpm.get_np(p_field))
    potexafmm = pfield_out[probe_i][vpm._n_param+1:vpm._n_param+3]
    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm, _ = vpm.conv(p_field, "ExaFMM")
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    # Stores particle values
    Us, convs = [], []
    for i in 1:n_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(part_Us, Us)
    push!(part_convs, convs)
    # Stores probe values
    Us, convs = [], []
    for i in n_particles+1:n_particles + n_probe_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(probe_Us, Us)
    push!(probe_convs, convs)

    if verbose
      println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$potexafmm")
        println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end
  end



  if exafmmcsda
    # ---------------- ExaFMM
    this_label = "CSDA ExaFMM"
    push!(method_labels, this_label)
    if verbose;println("\t***************** $this_label *****************");end;
    vpm.set_P2PTYPE(Int32(4))
    vpm.set_L2PTYPE(Int32(1))

    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "ExaFMM")
    _, _, _, _, pfield_out = vpm.U_ExaFMM(p_field._p_field, vpm.get_np(p_field))
    potexafmm = pfield_out[probe_i][vpm._n_param+1:vpm._n_param+3]
    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm, _ = vpm.conv(p_field, "ExaFMM")
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    # Stores particle values
    Us, convs = [], []
    for i in 1:n_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(part_Us, Us)
    push!(part_convs, convs)
    # Stores probe values
    Us, convs = [], []
    for i in n_particles+1:n_particles + n_probe_particles
      push!(Us, Mexafmm[i,4:6])
      push!(convs, Mexafmm[i,1:3])
    end
    push!(probe_Us, Us)
    push!(probe_convs, convs)

    if verbose
      println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$potexafmm")
        println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end
  end


    if exafmmregANAcsda
      # ---------------- ExaFMM singANA-CSDA
      this_label = "regANA-CSDA ExaFMM"
      push!(method_labels, this_label)
      if verbose;println("\t***************** $this_label *****************");end;
      vpm.set_P2PTYPE(Int32(3))
      vpm.set_L2PTYPE(Int32(1))

      if verbose; println("\tTest run..."); end;
      vpm.conv(p_field, "ExaFMM"; Uprobes=Uprobes)
      _, _, _, _, pfield_out = vpm.U_ExaFMM(p_field._p_field, vpm.get_np(p_field))
      potexafmm = pfield_out[probe_i][vpm._n_param+1:vpm._n_param+3]
      if verbose; println("\tNow for real..."); end;
      prev_t = time()
      Mexafmm, _ = vpm.conv(p_field, "ExaFMM"; Uprobes=Uprobes)
      end_t = time()
      texa = end_t-prev_t

      # Induced velocity
      Uexafmm = Mexafmm[probe_i,4:6]
      convexafmm = Mexafmm[probe_i,1:3]

      # Stores particle values
      Us, convs = [], []
      for i in 1:n_particles
        push!(Us, Mexafmm[i,4:6])
        push!(convs, Mexafmm[i,1:3])
      end
      push!(part_Us, Us)
      push!(part_convs, convs)
      # Stores probe values
      Us, convs = [], []
      for i in n_particles+1:n_particles + n_probe_particles
        push!(Us, Mexafmm[i,4:6])
        push!(convs, Mexafmm[i,1:3])
      end
      push!(probe_Us, Us)
      push!(probe_convs, convs)

      if verbose
        println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$potexafmm")
          println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
          println("\t\t\t$Uexafmm")
          println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
          println("\t\t\t$convexafmm")
        println("\tTime: $(end_t-prev_t) (s)")
      end
    end

  # --------------- ERROR
  if (exafmmcsda || exafmmfd) && (dirsing || dirblob)
    num = zeros(3)
    den = zeros(3)
    if error=="convective"
      for i in 1:n_particles
        num += (Mexafmm[probe_i,1:3] - Mdir[probe_i,1:3]).^2
        den += (Mdir[probe_i,1:3]).^2
      end
    else
      for i in 1:n_particles
        num += (Mexafmm[probe_i,4:6] - Mdir[probe_i,4:6]).^2
        den += (Mdir[probe_i,4:6]).^2
      end
    end
    rms = sqrt(num./den)
    if verbose && error=="convective"
      println("\n********** RMS error on the convective term field: $rms")
    elseif verbose
      println("\n********** RMS error on U field: $rms")
    end
  end


  # --------------- PLOT
  if plot_
    if verbose; println("Plotting..."); end;
    # PLOT
    fig = figure("particles")
    ax = gca(projection="3d")
    ## Particles
    x,y,z = [], [], []
    c_list = []
    for i in 1:n_particles+n_probe_particles
        this_pos = vpm.get_x(p_field, i)
        push!(x, this_pos[1])
        push!(y, this_pos[2])
        push!(z, this_pos[3])
        this_rgb = (i<=n_particles ? (vpm.get_gamma(p_field, i)-ql)./(qu-ql)
                                      : zeros(3))
        push!(c_list, this_rgb)
    end
    scatter3D(x, y, z, c=c_list)

    # Values plots
    dims_labels = ["x", "y", "z"]
    # dims_style = ["r", "g", "b", "c"]
    dims_style = ["", "", "", ""]
    # method_style = ["ko", "r*", "g^", "bx"]
    method_style = ["k.", "g.", "r.", "b.", "y.", "c.", "m."]

    ## Probes U
    fig = figure("probes", figsize=(7*3,5*1))
    for Ui in 1:3
      subplot(130 + Ui)
      for (i, label) in enumerate(method_labels)
        x = [p_i for p_i in 1:n_probe_particles]
        ref = [probe_Us[1][xi][Ui] for xi in x]
        y = [(probe_Us[i][xi][Ui]-ref[xi])/ref[xi] for xi in x]
        plot(x,y*100, method_style[i]*dims_style[Ui],
              label="$label")
        legend(loc="best")
        grid(true, color="0.8", linestyle="--")
        xlabel("Probe index")
        ylabel("U$(dims_labels[Ui]) error %")
        title("U$(dims_labels[Ui]) error relative to $(method_labels[1])")
      end
    end

    ## Particles U
    fig = figure("particlesvals", figsize=(7*3,5*2))
    crits = [[] for i in 1:size(method_labels)[1]]
    for Ui in 1:3
      subplot(230 + Ui)
      for (i, label) in enumerate(method_labels)
        x = [p_i for p_i in 1:n_particles]
        ref = [part_Us[1][xi][Ui] for xi in x]
        y = [(part_Us[i][xi][Ui]-ref[xi])/ref[xi] for xi in x]
        plot(x,y*100, method_style[i]*dims_style[Ui],
              label="$label")
        legend(loc="best")
        grid(true, color="0.8", linestyle="--")
        xlabel("Particle index")
        ylabel("U$(dims_labels[Ui]) error %")
        title("U$(dims_labels[Ui]) error relative to $(method_labels[1])")
        for (yi,yval) in enumerate(y)
          if yval>1 && false==(yi in crits[i]); push!(crits[i], yi); end;
        end
      end
    end
    println("------------ Criticals in U -----------")
    for (i,crits_meth) in enumerate(crits)
      sort!(crits_meth)
      for xi in crits_meth
        println("Part #$(xi)\n\t$(method_labels[1])=$(part_Us[1][xi])"*
                "\n\t$(method_labels[i])=$(part_Us[i][xi])"*
                "\n\tError=$((part_Us[1][xi]-part_Us[i][xi])./part_Us[1][xi])")
      end
    end
    ## Conv
    crits = [[] for i in 1:size(method_labels)[1]]
    for convi in 1:3
      subplot(230 + 3 + convi)
      for (i, label) in enumerate(method_labels)
        x = [p_i for p_i in 1:n_particles]
        ref = [part_convs[1][xi][convi] for xi in x]
        y = [(part_convs[i][xi][convi]-ref[xi])/ref[xi] for xi in x]
        plot(x,y*100, method_style[i]*dims_style[convi],
              label="$label")
        legend(loc="best")
        grid(true, color="0.8", linestyle="--")
        xlabel("Particle index")
        ylabel("conv$(dims_labels[convi]) error %")
        title("conv$(dims_labels[convi]) error relative to $(method_labels[1])")
        for (yi,yval) in enumerate(y)
          if yval>1 && false==(yi in crits[i]); push!(crits[i], yi); end;
        end
      end
    end
    println("------------ Criticals in conv -----------")
    for (i,crits_meth) in enumerate(crits)
      sort!(crits_meth)
      for xi in crits_meth
        println("Part #$(xi)\n\t$(method_labels[1])=$(part_convs[1][xi])"*
                "\n\t$(method_labels[i])=$(part_convs[i][xi])"*
                "\n\tError=$((part_convs[1][xi]-part_convs[i][xi])
                              ./part_convs[1][xi])")
      end
    end
  end

  if (exafmmcsda || exafmmfd) && (dirsing || dirblob)
    if error=="convective"
      err = (convdir - convexafmm)/convdir
    else
      err = (Udir - Uexafmm)/Udir
    end
  else
    tdir = nothing
    err = nothing
    texa = nothing
    rms = nothing
  end

  return tdir, texa, rms, err
end
