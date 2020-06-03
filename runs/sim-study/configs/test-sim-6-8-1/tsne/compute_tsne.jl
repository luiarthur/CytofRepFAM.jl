include("tsne_imports.jl")

get_simdat_path(path_to_output) = joinpath(splitdir(path_to_output)[1], "simdat.bson")
results_dir = "/scratchdata/alui2/cytof/results/repfam/test-sim-6-8-1"

for zind in 1:3
  for pmiss in (0.0, 0.6)
    # Get simulation name
    simname = "pmiss$(pmiss)-phi0.0-zind$(zind)"

    println("Processing: ", simname)

    # Get path to simuldate data.
    path = joinpath(results_dir, simname, "simdat.bson")

    # Fit TSNE jointly on all samples.
    tsne, sample_ind, Y, M = let
      use_complete_data = (pmiss == 0.0)
      compute_combined_tsne(path, use_complete_data=use_complete_data,
                            seed=0, verbose=2)
    end

    # Create data frame
    df = let
      num_markers = size(Y, 2)
      Ynames = [Symbol("Y$(j)") for j in 1:num_markers]
      Mnames = [Symbol("M$(j)") for j in 1:num_markers]
      colnames = [:tsne1, :tsne2, :sample_ind, Ynames..., Mnames...]
      DataFrame([tsne sample_ind Y M], colnames)
    end

    # Write output to file.
    csvdir = "out/csv"
    mkpath(csvdir)
    CSV.write("$(csvdir)/$(simname).csv", df)
  end
end



