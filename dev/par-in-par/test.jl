using Distributed

addprocs(12)

@everywhere include("myfunc.jl")

_ = myfunc(1)
@time w = myfunc(1)  # 6.6s

_ = mypfunc(2)
@time x = mypfunc(2)  # 9.6s

# NOTE: This version is SLOW!
_ = pmap(mypfunc, 1:2)
@time y = pmap(mypfunc, 1:2)  # 16.1s

# NOTE: This seems like the way to go.
_ = asyncmap(mypfunc, 1:2)
@time z = asyncmap(mypfunc, 1:2)  # 11.7s

@time @sync begin
  for i in 1:2
    @async mypfunc(i)
  end
  println("DONE")
end  # 11.9s

rmprocs(workers())

# times = [@elapsed myfunc(0) for _ in 1:20]
