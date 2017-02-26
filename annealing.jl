cd("/home/oskar/Documents/BU/Research/Annealing")

module annealing
  push!(LOAD_PATH, pwd())
  import Graphslib
  Gl = Graphslib

  ∧(A::Bool, B::Bool) = A && B
  ∨(A::Bool, B::Bool) = A || B

  function idxtopos(i::Int, ds::Array{Int})
    x = round(Int, i/(ds[2]), RoundUp)
    y = i - ds[2] * (x-1)
    return [x,y]
  end

  function postoidx(pos::Array{Int}, ds::Array{Int})
    return (pos[1]-1) * ds[2] + pos[2]
  end

  # 2-bit gates
  function swap(i::Int, swap1::Int=1, swap2::Int=2, reverse::Bool=false)
    a = bits(i)
    a = a[1:end-swap2] * string(a[end-swap1+1]) * string(a[end-swap2+1]) * a[end-swap1+2:end]
    parse(Int,a,2)
  end
  function iden(i::Int, reverse::Bool=false) i end

  # 3-bit gates
  function toff(i::Int, reverse::Bool=false)
    bits(i)[end-1:end] == "11" && return i $ 4
    return i
  end
  function idid(i::Int, reverse::Bool=false)  i end
  function swid(i::Int, reverse::Bool=false)  swap(i, 2, 3) end
  function idsw(i::Int, reverse::Bool=false)  swap(i, 1, 2) end
  function swsw(i::Int, reverse::Bool=false)
    reverse ? return swap(swap(i,1,2),2,3) : return swap(swap(i,2,3),1,2)
  end

  function random_gatemap(ds::Array{Int}, ps::Array)
    ps /= sum(ps)
    ngates = prod(ds)
    gm = zeros(Int, ngates)
    for i in 1:length(gm)
      check = 0
      random_number = rand()
      for j in 1:length(ps)
        check += ps[j]
        if random_number < check
          gm[i] = j
          break
        end
      end
    end
    return gm
  end

  function vertex_lattice_gategraph(ds::Array{Int})
    ngates = prod(ds)
    g = Gl.newgraph(ngates)
    g.vs["state"] = -1 * ones(Int, ngates)

    for i in 1:ngates
      pos = idxtopos(i,ds)
      g.vs[i]["x"] = pos[1]
      g.vs[i]["y"] = pos[2]
    end

    edge_index = 0
    for vertex1_index in 1: ngates - ds[2]
      pos = idxtopos(vertex1_index,ds)
      vertex2_index = postoidx(pos + [1,0], ds)
      Gl.add_edges!(g, vertex1_index, vertex2_index)

      edge_index += 1

      if mod(pos[1], 2) ≡ 1
        g.es[edge_index]["bits"] = [[3,2],[2,1]]
      else
        g.es[edge_index]["bits"] = [[1],[3]]
      end

      dy = (-1)^(mod(pos[1],2))
      y = pos[2] + dy

      if y ≡ ds[2]+1
        y=1
      elseif y ≡ 0
        y = ds[2]
      end
      vertex2_index = postoidx( [pos[1]+1, y], ds)
      Gl.add_edges!(g, vertex1_index, vertex2_index)

      edge_index += 1

      if mod(pos[1], 2) ≡ 1
        g.es[edge_index]["bits"] = [[1],[3]]
      else
        g.es[edge_index]["bits"] = [[3,2],[2,1]]
      end
    end
    return g
  end

  function vertex_lattice(ds::Array{Int}, gates::Array, nbits, gm::Array{Int}=Array{Int}(0), st::Array{Int}=Array{Int}(0), ps::Array{AbstractFloat}=Array{AbstractFloat}(0))

    if typeof(nbits) ≡ Int
      nbits = nbits * ones(Int, length(gates))
    elseif length(nbits) ≠ length(gates)
      error("length(nbits) ≠ length(gates)")
    elseif !isempty(ps) ∧ length(ps) ≠ length(gates)
      error("length(ps) ≠ length(gates)")
    elseif !isempty(ps) ∧ !isempty(gm)
      error("both ps and gm specified")
    end

    if !isempty(ps)
      gm = random_gatemap(ds, ps)
    elseif isempty(gm)
      ps = rand(length(gates))
      ps /= sum(ps)
      gm = random_gatemap(ds, ps)
    end
    g = vertex_lattice_gategraph(ds)

    g.vs["type"] = gm
    for i in 1:length(g.vs)
      g.vs[i]["gate"] = gates[g.vs[i]["type"]]
      g.vs[i]["nbits"] = nbits[g.vs[i]["type"]]
    end
    isempty(st) || g.vs["state"] = st
    # g.vs["tensor"] = init_tensors(g)
    return g
  end

  function compute_energy(gate1::Function, gate2::Function, bits1::Array{Int}, bits2::Array{Int}, state1::Int, state2::Int)
    energy = 0
    for i in 1:length(bits1)
      energy += bits(gate1(state1))[end-bits1[i]+1] !== bits(gate2(state2))[end-bits2[i]+1]
    end
    energy
  end

  function solve_gategraph!(g::Gl.Graph)

    for vertex in g.vs
      vertex["state"] ≠ -1 && continue
      neighbors_index = Gl.neighbors(g, vertex, true, false)
      isempty(neighbors_index) && continue
      -1 in g.vs[neighbors_index]["state"] && continue
      gatebits = Array{String}(vertex["nbits"])
      for neighbor in g.vs[neighbors_index]
        edge_index = Gl.get_eid(g, neighbor.index, vertex.index)
        gatebits[g.es[edge_index]["bits"][2]] = split(bits(neighbor["gate"](neighbor["state"]))[end+1-g.es[edge_index]["bits"][1]], "")
      end
      vertex["state"] = parse(Int, join(gatebits[end:-1:1]),2)
    end
    for vertex in g.vs[end:-1:1]
      vertex["state"] ≠ -1 && continue
      neighbors_index = Gl.neighbors(g, vertex, false, true)
      isempty(neighbors_index) && continue
      -1 in g.vs[neighbors_index]["state"] && continue
      gatebits = Array{String}(vertex["nbits"])
      for neighbor in g.vs[neighbors_index]
        edge_index = Gl.get_eid(g, neighbor.index, vertex.index)
        gatebits[g.es[edge_index]["bits"][1]] = split(bits(neighbor["gate"](neighbor["state"]))[end+1-g.es[edge_index]["bits"][2]], "")
      end
      vertex["state"] = parse(Int,join(gatebits[end:-1:1]),2)
    end
    return g.vs["state"]
  end

  function check_graph_solution(g::Gl.Graph)
    errnumber = 0
    for edge in g.es
      if g.vs[edge.source]["state"] ≡ -1 || g.vs[edge.target]["state"] ≡ -1 continue end
      if bits(g.vs[edge.target]["state"])[end+1-edge["bits"][2]] ≠ bits(g.vs[edge.source]["gate"](g.vs[edge.source]["state"]))[end+1-edge["bits"][1]]
        println("Edge $(edge.index) from gate $(edge.source) to gate $(edge.target) doesn't match.")
        errnumber += 1
      end
    end
    errnumber ≠ 0 && println("$errnumber edges don't match in the graph.")
  end
end
#################################################################################################################################################################

function init_annealing!(g::Gl.Graph, ds::Array{Int}, fixed_left::Array{Int}, fixed_right::Array{Int})
  g.vs["state"] = rand(1:8, length(g.vs["state"]))
  g.vs["fixed"] = false
  g.vs[fixed_left]["fixed"] = true
  g.vs[end - ds[2] + fixed_right]["fixed"] = true
  return
end


function anneal!(g::Gl.Graph, mcsteps::Int=1000, βₘᵢₙ::AbstractFloat=0., βₘₐₓ::AbstractFloat=4., majority_steps::Int=100, majority_cutoff::AbstractFloat=0.9)
  Δβ = (βₘₐₓ - βₘᵢₙ) / mcsteps
  States = Array{Int}(length(g.vs), majority_steps)
  for j in 1:majority_steps
    β = βₘᵢₙ
    for i in 1:mcsteps
      montecarlo(g, β)
      β += Δβ
    end
    States[:,j] = g.vs["state"]
  end
  #majority(g, States, majority_cutoff, majority_steps)
end

function montecarlo!(g::Gl.Graph, β::AbstractFloat)
  for i in 1:length(g.vs)
    Vert = g.vs[rand(1:length(g.vs))]
    Vert["fixed"] && continue
    neighbor_vertices_index = Gl.neighbors(g, Vert)
    ΔE = 0
    for neighbor_Vert_index in neighbor_vertices_index
      neighbor_Vert = g.vs[neighbor_Vert_index]
      edge = g.es[Gl.get_eid(g, Vert.index, neighbor_Vert.index)[1]]
      if edge.source == Vert.index
        ΔE -= an.compute_energy(Vert["gate"], neighbor_Vert["gate"], edge["bits"][1], edge["bits"][2], Vert["state"], neighbor_Vert["state"])
      else
        ΔE -= an.compute_energy(Vert["gate"], neighbor_Vert["gate"], edge["bits"][2], edge["bits"][1], Vert["state"], neighbor_Vert["state"])
      end
    end
      newstate = rand(1:8)

    for neighbor_Vert_index in neighbor_vertices_index
      neighbor_Vert = g.vs[neighbor_Vert_index]
      edge = g.es[Gl.get_eid(g, Vert.index, neighbor_Vert.index)[1]]
      if edge.source == Vert.index
        ΔE += an.compute_energy(Vert["gate"], neighbor_Vert["gate"], edge["bits"][1], edge["bits"][2], newstate, neighbor_Vert["state"])
      else
        ΔE += an.compute_energy(Vert["gate"], neighbor_Vert["gate"], edge["bits"][2], edge["bits"][1], newstate, neighbor_Vert["state"])
      end
      rand() < exp(- β * ΔE) && Vert["state"] = newstate
    end
  end
end

function majority(g::Gl.Graph, States::Array{Int}, majority_cutoff::AbstractFloat, majority_steps::Int)
  for i in 1:length(g.vs)
    state_counts = counts(States[i,:], 1:8) / majority_steps
    for (j, state_count) in enumerate(state_counts)
      if state_count > majority_cutoff
        g.vs[i]["state"] = j
        g.vs[i]["fixed"] = true
      end
    end
  end
end


import StatsBase.counts



Gl = Graphslib
using annealing
an = annealing
supertype(an.idid)
gates = [an.idid]
ds = [40,20]
32*4
@time anneal!(g, 2^10, 0., 8., 1, 0.9)
g = an.vertex_lattice(ds,gates, 3)
g.vs[fixed_left]["fixed"] = true
fixed_right = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
fixed_left = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
init_annealing!(g, ds, fixed_left, fixed_right)
g.vs["fixed"]
g.vs["fixed"] = false
g.vs[1].index
g.vs["state"] = Array{Int}(g.vs["state"])
g.vs["state"]
print(g.vs["fixed"])
Gl.get_eid(g,1,6)
g = vertex_lattice_gategraph([10,5])
g.vs["type"] = gm
ps =rand(5)
gates[g.vs[1]["type"]]
ds = [10,5]
a = [1]
gm = random_gatemap(ds,ps)
g.vs[1]["type"]
length(g.vs) == length(gm)
end
Pkg.add("GLVisualize.jl"))
