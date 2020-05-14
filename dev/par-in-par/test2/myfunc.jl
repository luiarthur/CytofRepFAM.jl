function myfunc(n)
  println("inner pid: $(getpid())")
  out = 0.0
  x = @elapsed for i in 1:Int(2e5)
    _ = sum(ones(10000) .^2)
  end

  return x
end

function mypfunc(n)
  # NOTE: The redirects cause problems.
  # One fix would be to explicitly, always write to file in src code
  # redirect_stdout_to_file("$(n).txt") do
    println("Start $(n)! pid: $(getpid())")
    x = pmap(myfunc, collect(1:3) .+ n)
    println("Finished $(n)!")
    return x
  # end
end

 
function redirect_stdout_to_file(f::Function, path::String)
  open(path, "w") do io
    redirect_stdout(io) do
      f()
    end
  end
end

