function myfunc(n)
  # println(n)
  out = 0.0
  # x = @elapsed for i in 1:1000
  #   n = 2000000
  #   out += sum(randn(n) .^ 2) / n
  # end
  x = @elapsed for i in 1:Int(2e5)
    _ = sum(ones(10000) .^2)
  end

  return x
end

function mypfunc(n)
  println("Start $(n)!")
  x = pmap(myfunc, collect(1:6) .+ n)
  println("Finished $(n)!")
  return x
end
