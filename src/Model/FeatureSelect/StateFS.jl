"""
Parameter state in Gibbs sampler for feature allocation model
with feature selection. This augments `State` with parameters:
`r`, `W_star`, and `p`-or-`omega`.

Usage example:
==============

```julia
julia> s = StateFS{Float64}()
julia> s.r = rand(Bool, 3, 10)
```

NOTE:
=====
The update order should be:
Z -> v -> alpha -> p -> r -> lam -> W* -> gamma -> eta -> delta -> sig^2 -> y*

Particularly, `lam` must be updated before W* and preferrably after p & r.
And `gamma` should be updated after Z and before delta.
"""
mutable struct StateFS{F <: AbstractFloat}
  r::Matrix{Bool}  # dim: IxK
  W_star::Matrix{F}  # dim: IxK
  omega::Vector{F}  # dim: number of covariates (`length(x_i)`) + 1
  theta::State{F}  # State

  # Primary constructor returns an uninitialized stateFS.
  StateFS{F}() where {F <: AbstractFloat} = new()
end


function StateFS{F}(theta::State{F}, dfs::DataFS;
                    omega_prior=nothing,
                    eps_r::Float64=0.0, verbose=1) where {F <: AbstractFloat}
  sfs = StateFS{F}()
  I, K = size(theta.W)

  # Set theta
  sfs.theta = theta

  # Set W*
  sfs.W_star = theta.W * 5

  # Set R
  if eps_r > 0.0
    sfs.r = sfs.theta.W .> eps_r
  else
    sfs.r = ones(Bool, I, K)
  end

  # Set omega
  if omega_prior == nothing
    sfs.omega = zeros(dfs.P)
  else
    sfs.omega = rand(omega_prior, dfs.P)
  end

  if verbose > 0
    println("Initial sfs.r: $(sfs.r)")
    println("Initial sfs:W_star $(sfs.W_star)")
    println("Initial sfs.theta.W: $(sfs.theta.W)")
    println("Initial sfs.theta.Z: $(sfs.theta.Z)")
    println("Inital sfs.omega: $(sfs.omega)")
  end

  return sfs
end
