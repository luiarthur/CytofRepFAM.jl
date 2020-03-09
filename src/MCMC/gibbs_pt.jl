"""
TODO...
"""
function gibbs_pt(init::Vector{T}, args::Vector{Dict{Symbol, Any}},
                  update::Function;
                  monitors::Vector{Vector{Symbol}}=deepcopy(monitor_default),
                  thins::Vector{Int}=deepcopy(thin_default),
                  nmcmc::Int64=1000, nburn::Int=0,
                  printFreq::Int=0,
                  save_all_states::Bool=false,
                  printlnAfterMsg::Bool=true) where T

  @assert printFreq >= -1
  if printFreq == 0
    numPrints = 10
    printFreq = Int(ceil((nburn + nmcmc) / numPrints))
  end

  states = deepcopy(init)
  num_chains = length(states)

  # Checking number of monitors.
  numMonitors = length(monitors)
  println("Number of monitors: $(numMonitors)"); flush(stdout)
  @assert numMonitors == length(thins)

  # Check monitor
  if numMonitors == 0
    println("Using default monitor."); flush(stdout)
    fnames = [fname for fname in fieldnames(typeof(init))]
    append!(monitors, [fnames])
    append!(thins, 1)
    numMonitors = 1
  end

  # Number of Samples for each Monitor
  numSamps = [div(nmcmc, thins[i]) for i in 1:numMonitors]

  if printFreq > 0
    println("Preallocating memory..."); flush(stdout)
  end

  # Create object to return
  if save_all_states
    @time out = [fill([deepcopyFields(s, monitors[i]) for s in states],
                      numSamps[i]) for i in 1:numMonitors]
  else
    @time out = [fill(deepcopyFields(states[1], monitors[i]), numSamps[i])
                 for i in 1:numMonitors]
  end

  function printMsg(i::Int)
    if (printFreq > 0) && (i % printFreq == 0) && (i > 1)
      if :ll in keys(args[1])
        loglikeMsg = "-- loglike: $(last(args[1][:ll]))"
      else
        loglikeMsg = ""
      end

      print("$(showtime()) -- $(i)/$(nburn+nmcmc) $loglikeMsg")

      if printlnAfterMsg
        println()
      end

      flush(stdout)
    end
  end

  # burn in
  for i in 1:nburn
    states, args = update(states, args, i)
    printMsg(i)
  end

  counters = zeros(Int, numMonitors)
  # Gibbs loop
  for i in 1:nmcmc
    states, args = update(states, args, i + nburn)
    printMsg(i + nburn)

    for j in 1:numMonitors
      if i % thins[j] == 0
        if save_all_states
          substate = [deepcopyFields(s, monitors[j]) for s in states]
        else
          substate = deepcopyFields(states[1], monitors[j])
        end
        counters[j] += 1
        out[j][counters[j]] = substate
      end
    end
  end

  return (out, states, args)
end
