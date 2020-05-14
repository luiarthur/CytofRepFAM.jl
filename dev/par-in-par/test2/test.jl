# using Distributed

# addprocs(8)
# @everywhere include("myfunc.jl")
# 
# 
# # MAIN
# let
#   buffer = Task[]
#   delay = 1.0
#   blocks = 2
# 
#   for i in 1:6
#     t = @async mypfunc(i)
#     append!(buffer, [t])
# 
#     if i % blocks == 0
#       foreach(wait, buffer)
#       buffer = Task[]
#     end
#   end
# end
# 
# rmprocs(workers())

# MAIN
using Distributed
addprocs(8)
@everywhere include("myfunc.jl")

@time asyncmap(mypfunc, 1:12, ntasks=2)
# rmprocs(workers())
