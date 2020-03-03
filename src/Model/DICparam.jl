@namedargs mutable struct DICparam
  p::Vector{Matrix{Float64}}
  mu::Vector{Matrix{Float64}}
  sig::Vector{Vector{Float64}}
  y::Vector{Matrix{Float64}}
end


function computeLoglikeDIC(d::Data, param::DICparam)
  ll = 0.0

  for i in 1:d.I
    for j in 1:d.J
      for n in 1:d.N[i]
        y_inj_is_missing = (d.m[i][n, j] == 1)

        # NOTE: Refer to `../compute_loglike.jl` for reasoning.
        if y_inj_is_missing

          # Compute p(m_inj | y_inj, theta) term.
          ll += log(param.p[i][n, j])
        end

        # Compute p(y_inj | theta) term.
        ll += logpdf(Normal(param.mu[i][n, j], param.sig[i][n]),
                     param.y[i][n, j])
      end
    end
  end

  return ll
end

function updateParams(d::MCMC.DICstream{DICparam}, param::DICparam)
  d.paramSum.p += param.p
  d.paramSum.mu += param.mu
  d.paramSum.sig += param.sig
  d.paramSum.y += param.y

  return
end

function paramMeanCompute(d::MCMC.DICstream{DICparam})::DICparam
  return DICparam(d.paramSum.p / d.counter,
                  d.paramSum.mu / d.counter,
                  d.paramSum.sig / d.counter,
                  d.paramSum.y / d.counter)
end
