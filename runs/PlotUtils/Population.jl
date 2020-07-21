mutable struct Population
  all::Dict{String, Int}
  size::Int
end

Population() = Population(Dict(), 0)
Base.keys(p::Population) = keys(p.all)

function binarize(p::Population, subpop)
  subpop_char_list = [string(Int(s)) for s in subpop]
  return join(subpop_char_list, "")
end

function label(p::Population, subpop)
  b = binarize(p, subpop)
  if !(b in keys(p))
    p.size += 1
    p.all[b] = p.size
  end
  return p.all[b]
end

# Example usage:
# p = Population()
# label(p, [1,0,1])  # 1
# label(p, [1,0,1])  # 1
# label(p, [1,1,1])  # 2
# label(p, [1,1,1,1]) # 3
