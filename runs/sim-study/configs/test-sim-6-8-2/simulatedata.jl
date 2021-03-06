include(joinpath(@__DIR__, "../../rand_skewnormal.jl"))

function simulatedata(; Z, W, N::Vector{Int},
                      mus::Dict{Int, Vector{Float64}},
                      sig2=Vector{Float64},
                      skew::Float64,
                      sortLambda=false, seed=nothing,
                      propmissingscale=0.3,
                      expressed_not_missing=false,
                      beta=[-9.2, -2.3],  # linear missing mechanism
                      eps_mus_dist=Uniform(-.3, .3))
  if seed != nothing
    Random.seed!(seed)
  end

  J, K = size(Z)
  I = length(N)

  # Create dict of L
  L = Dict(0 => length(mus[0]),
           1 => length(mus[1]))

  @assert size(W) == (I, K)
  @assert length(sig2) == I
  @assert all(sig2 .> 0)
  @assert all(length(mus[key]) == L[key] for key in keys(L))
  @assert all(N .> 0)

  # Simulate cell phenotypes (lambda)
  lam = [begin
           rand(Categorical(W[i, :]), N[i])
         end for i in 1:I]

  if sortLambda
    lam = [sort(lami) for lami in lam]
  end

  # Generate y
  y = [zeros(Ni, J) for Ni in N]
  y_complete = deepcopy(y)

  z(i::Int, n::Int, j::Int) = Z[j, lam[i][n]]

  eps_mus = let
    if eps_mus_dist != nothing
      Dict(0 => rand(eps_mus_dist, I, J),
           1 => rand(eps_mus_dist, I, J))
    else
      Dict(0 => zeros(I, J),
           1 => zeros(I, J))
    end
  end

  function mu(i::Int, n::Int, j::Int)
    @assert L[0] == L[1] == 1
    z_inj = z(i, n, j)
    return mus[z_inj][1] + eps_mus[z_inj][i, j]
  end

  for i in 1:I
    sig_i = sqrt(sig2[i])
    for j in 1:J
      for n in 1:N[i]
        # y_complete[i][n, j] = rand(Normal(mu(i, n, j), sig_i))
        y_complete[i][n, j] = rand_skewnormal(mu(i, n, j), sig_i, skew)
      end
      # Set some y to be missing
      # Using a linear missing mechanism to rate missingness
      p_miss = [CytofRepFAM.Model.prob_miss(y_complete[i][n, j], beta)
                for n in 1:N[i]]

      # If I want to truly expressed cells to not be missing:
      if expressed_not_missing
        # If Z_inj == 1, then p_miss = 0
        # Otherwise, p_miss = p_miss
        p_miss .*= [1 - z(i, n, j) for n in 1:N[i]]
      end

      # If this is not used, then many positive markers may be missing too.
      prop_not_expressed = sum(W[i,:] .* (1 .- Z[j,:]))
      # If propmissingscale is 1, about all the non-expressed observations
      # will be missing. Best to set propmissingscale to about 0.3.
      prop_missing = propmissingscale * prop_not_expressed

      num_missing = Int(round(N[i] * prop_missing))
      idx_missing = Distributions.wsample(1:N[i], p_miss, num_missing,
                                          replace=false)
      y[i][:, j] .= y_complete[i][:, j] .+ 0

      if propmissingscale > 0
        y[i][idx_missing, j] .= NaN
      end
    end
  end

  return Dict(:Z => Z, :N => N, :L => L, :mus => mus, :W => W,
              :seed => seed, :lam => lam, :sig2 => sig2, :y => y,
              :y_complete => y_complete, :beta => beta,
              :eps_mus => eps_mus, :skew => skew)
end
