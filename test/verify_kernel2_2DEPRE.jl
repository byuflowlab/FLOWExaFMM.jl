VPM_PATH = "/home/user/Dropbox/FLOWResearch/MyCodes/MyVPM/src/"
# VPM_PATH = "MyVPM/"
include(VPM_PATH*"MyVPM.jl")
vpm = MyVPM

using PyPlot
pygui(false)

function run(n_particles; probe_i=1, verbose=true, plot=false,
                    P::Int32=vpm.exafmm_P, NCRIT::Int32=vpm.exafmm_NCRIT,
                    THETA::Float64=vpm.exafmm_THETA, fmmverbose=false,
                    error="velocity", sigmablob=0.1)
  vpm.set_P(P)
  vpm.set_NCRIT(NCRIT)
  vpm.set_THETA(THETA)
  vpm.set_VERBOSE(fmmverbose)

  len = 50                    # Length of box edge
  C = [100, -50, 25]          # Center of box
  l = C - len*[0.5,0.5,0.5]   # Lower bound for spawning particles
  u = C + len*[0.5,0.5,0.5]   # Upper bound for spawning particles

  qlen = 1                    # Dispersion of charge
  qC = [1, 1, 1]              # Center of charge
  ql = qC - qlen*[0.5,0.5,0.5]# Lower bound of charges
  qu = qC + qlen*[0.5,0.5,0.5]# Upper bound of charges

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
  p_field = vpm.ParticleField(n_particles, Uinf, nothing, "directblob")
  for i in 1:n_particles
    vpm.addparticle(p_field, vcat(Q[i], q[i], sigmablob) )
  end

  # Probe
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


    # ---------------- DEPRICATED ExaFMM
    if verbose;println("\t***************** DEPRICATED ExaFMM *****************");end;
    vpm.set_CSDA(false)
    vpm.set_REGULARIZED(true)
    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "DEPRIExaFMM")

    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm = vpm.conv(p_field, "DEPRIExaFMM")
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    if verbose
      println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end


  # ------------ MyVPM
  if verbose;println("\t***************** DirectBlob ***************************");end;

  if verbose; println("\tTest run..."); end;
  vpm.conv(p_field, "direct")

  if verbose; println("\tNow for real..."); end;
  prev_t = time()
  Mdir = vpm.conv(p_field, "direct")
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

  if verbose
    println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
    println("\t\t\t$potdir")
      println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$Udir")
      println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$convdir")
    println("\tTime: $(end_t-prev_t) (s)")
  end



    # ---------------- ExaFMM
    if verbose;println("\t***************** ExaFMM *****************************");end;
    vpm.set_CSDA(true)
    vpm.set_REGULARIZED(false)

    if verbose; println("\tTest run..."); end;
    vpm.conv(p_field, "ExaFMM")
    _, _, _, _, pfield_out = vpm.U_ExaFMM(p_field._p_field, vpm.get_np(p_field))
    potexafmm = pfield_out[probe_i][vpm._n_param+1:vpm._n_param+3]
    if verbose; println("\tNow for real..."); end;
    prev_t = time()
    Mexafmm = vpm.conv(p_field, "ExaFMM")
    end_t = time()
    texa = end_t-prev_t

    # Induced velocity
    Uexafmm = Mexafmm[probe_i,4:6]
    convexafmm = Mexafmm[probe_i,1:3]

    if verbose
      println("\t\tPotential at X=$(vpm.get_x(p_field,probe_i))")
      println("\t\t\t$potexafmm")
        println("\t\tVelocity induced at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$Uexafmm")
        println("\t\tConvective term at X=$(vpm.get_x(p_field,probe_i))")
        println("\t\t\t$convexafmm")
      println("\tTime: $(end_t-prev_t) (s)")
    end


  # --------------- ERROR
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



  # --------------- PLOT
  if plot
    if verbose; println("Plotting..."); end;
    # PLOT
    fig = figure()
    ax = gca(projection="3d")
    ## Particles
    x,y,z = [], [], []
    c_list = []
    for i in 1:n_particles
        this_pos = vpm.get_x(p_field, i)
        push!(x, this_pos[1])
        push!(y, this_pos[2])
        push!(z, this_pos[3])
    end
    scatter3D(x, y, z)
  end

  if error=="convective"
    err = (convdir - convexafmm)/convdir
  else
    err = (Udir - Uexafmm)/Udir
  end

  return tdir, texa, rms, err
end
