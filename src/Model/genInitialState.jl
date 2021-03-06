function genInitialState(c::Constants, d::Data;
                         verbose=1, sb_ibp=true,
                         allow_repeated_Z_columns=true)
  if verbose > 0
    println("Doing random initialization ....")
  end

  J = d.J
  K = c.K
  L = c.L
  I = d.I
  N = d.N

  vec_y = vcat(vec.(d.y)...)
  y_neg = filter(y_inj -> !isnan(y_inj) && y_inj < 0, vec_y)

  y_imputed = begin
    local out = [zeros(Float64, N[i], J) for i in 1:I]
    grid_size = 30
    y_grid = collect(range(-7, 0, length=grid_size))
    miss_probs = [[prob_miss(yg, c.beta[:, i]) for yg in y_grid] for i in 1:I]
    for i in 1:I
      for n in 1:N[i]
        for j in 1:J
          if isnan(d.y[i][n, j])
            # out[i][n, j] = rand(Uniform(y_lower, y_upper))
            out[i][n, j] = wsample(y_grid, miss_probs[i])
          else
            out[i][n, j] = d.y[i][n, j]
          end
          @assert !isnan(out[i][n, j])
        end
      end
    end

    out
  end

  alpha = mean(c.alpha_prior)
  v = [mean(Beta(alpha, 1)) for k in 1:K]
  b = sb_ibp ? cumprod(v) : v
  init_Z() = [Bool(rand(Bernoulli(b[k]))) for j in 1:J, k in 1:K]
  Z = init_Z()
  if !allow_repeated_Z_columns
    num_unique_cols(Z) = size(unique(Z, dims=2), 2)

    while num_unique_cols(Z) < K
      Z = init_Z()
    end
  end
  delta = Dict(false => rand(c.delta_prior[0], L[0]),
               true  => rand(c.delta_prior[1], L[1]))
  sig2 = [rand(c.sig2_prior) for i in 1:I]
  W = Matrix{Float64}(hcat([rand(c.W_prior) for i in 1:I]...)')
  lam = [ Int8.(rand(Categorical(W[i, :]), N[i])) for i in 1:I ]
  eta = begin
    function gen(z)
      arrMatTo3dArr([ rand(c.eta_prior[z]) for i in 1:I, j in 1:J ])
    end
    Dict([Bool(z) => gen(z) for z in 0:1])
  end
  gam = [zeros(Int8, N[i], J) for i in 1:I]
  for i in 1:I
    for j in 1:J
      for n in 1:N[i]
        z_lin = Z[j, lam[i][n]]
        gam[i][n, j] = rand(Categorical(eta[z_lin][i, j, :]))
      end
    end
  end

  eps = mean.(c.eps_prior)

  # Create State object
  state = State{Float64}()

  # Initialize State
  state.Z = Z
  state.delta = delta
  state.alpha = alpha
  state.v = v
  state.W = W
  state.sig2 = sig2
  state.eps = eps
  state.eta = eta
  state.lam = lam
  state.gam = gam
  state.y_imputed = y_imputed

  return state
end

include("SmartInit.jl")
