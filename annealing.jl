"""
  Random vertex lattice module
  First instance February 2017 by Oskar Pfeffer

"""
module annealing
  push!(LOAD_PATH, pwd())
  using Graphslib
  using JLD
  import StatsBase.counts
  export idxtopos, postoidx, swap, iden, idid, idsw, swid, toff, swsw, random_gatemap, edgebitsmap, vertex_lattice_gategraph, vertex_lattice, compute_energydiff, solve_gategraph!, solve_gategraph, check_graph_solution, init_annealing!, anneal!, montecarlo!, majority, totalunfits, unique_solution, energymeasure, statecounts, get_uniquesolutions

  # Faster and nicer definition for AND and OR
  ∧(A::Bool, B::Bool) = A && B
  ∨(A::Bool, B::Bool) = A || B

  function idxtovec(i::Int,l::Array{Int})
    """ Transforms scalar index to vector index. """
    ( i > prod(l)) && (return nothing)
    ( i < 1 ) && (return nothing)
    a = zeros(Int, (length(l)+1))
    idxtovec = zeros(Int, length(l))
    a[1] = 1
    a[2:end] = l
    for j in length(l):-1:1
      idxtovec[j] = fld((i-1),prod(a[1:j])) % a[j+1] + 1
    end
    return idxtovec
  end

  function vectoidx(r::Array{Int},l::Array{Int})
    """ Transforms vector index to scalar index. """
    ( length(l) != length(r) ) && (return nothing)
    vectoidx = 0
    a = zeros(Int, length(l)+1)
    a[1] = 1
    a[2:end] = l
    for i in length(l):-1:1
      ( r[i] > l[i] ) && (return nothing)
      ( r[i] < 1 ) && (return nothing)
      vectoidx += prod(a[1:i])*(r[i]-1)
    end
    vectoidx + 1
  end

  function idxtopos(i::Int,ds::Array{Int})
    """ Index to position vector. Index along y first. """
    return (idxtovec(i,ds[end:-1:1]))[end:-1:1]
  end

  function postoidx(pos::Array{Int}, ds::Array{Int})
    """ Position vector to index. Index along y first. """
    return vectoidx(Array(pos)[end:-1:1],ds[end:-1:1])
  end

  # 2-bit gates
  function swap(i::Int, swap1::Int=1, swap2::Int=2, reverse::Bool=false)
    a = bin(i,swap2)[end-swap2+1:end]
    a = string(a[end-swap1+1]) * string(a[end-swap2+1]) * string(a[end-swap1+2:end])
    parse(Int,a,2)
  end
  function iden(i::Int, reverse::Bool=false) i end

  # 3-bit gates
  function toff(i::Int, reverse::Bool=false)
    bits(i)[end-1:end] == "11" && return i $ 4
    return i
  end
  function idid(i::Int, reverse::Bool=false)  i end
  function swid(i::Int, reverse::Bool=false)
    i==2 && (return 4)
    i==4 && (return 2)
    i==3 && (return 5)
    i==5 && (return 3)
    return i
    #return swap(i, 2, 3)
  end
  function idsw(i::Int, reverse::Bool=false)
    i==1 && (return 2)
    i==2 && (return 1)
    i==5 && (return 6)
    i==6 && (return 5)
    return i
    #return swap(i, 1, 2)
  end
  function swsw(i::Int, reverse::Bool=false)
    if reverse
      i==1 && (return 4); i==2 && (return 1); i==3 && (return 5)
      i==4 && (return 2); i==5 && (return 6); i==6 && (return 3)
      return i
    else
      i==1 && (return 2); i==2 && (return 4); i==3 && (return 6)
      i==4 && (return 1); i==5 && (return 3); i==6 && (return 5)
      return i
    end
    #reverse ? swap(swap(i,1,2),2,3) : swap(swap(i,2,3),1,2)
  end

  function random_gatemap(ds::Array{Int}, ps::Array)
    """
    Returns a random gatemap (Integer from 1 to length(ps)) for a rectangular
    n-dimensional lattice with given probabilities ps for each gate
    """
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
  """Returns the map of the inputbits and/or outputbits of the edges of the graph g."""
  function edgebitsmap(g::Graph; input::Bool=true, output::Bool=true)
    edgebitsoutmap = Vector{Vector{Int}}(length(g.es))
    edgebitsinmap = Vector{Vector{Int}}(length(g.es))
    for edge in g.es
      edgebitsoutmap[edge.index] = edge["bits"][1]
    end
    for edge in g.es
      edgebitsinmap[edge.index] = edge["bits"][2]
    end
    (input && output) && return (edgebitsinmap, edgebitsoutmap)
    input && return edgebitsinmap
    output && return edgebitsoutmap
  end


  function vertex_lattice_gategraph(ds::Array{Int})
    """
    Returns an initialized graph with connected vertices for a 2-dimensional
    vertex lattice.
    """
    ngates = prod(ds); g = newgraph(ngates)
    g.vs["state"] = -1 * ones(Int, ngates)
    for i in 1:ngates
      pos = idxtopos(i,ds)
      g.vs[i]["x"] = pos[1]; g.vs[i]["y"] = pos[2]
    end
    edge_index = 0
    for vertex1_index in 1: ngates - ds[2]
      pos = idxtopos(vertex1_index,ds)
      vertex2_index = postoidx(pos + [1,0], ds)
      add_edges!(g, vertex1_index, vertex2_index)
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
      add_edges!(g, vertex1_index, vertex2_index)
      edge_index += 1
      if mod(pos[1], 2) ≡ 1
        g.es[edge_index]["bits"] = [[1],[3]]
      else
        g.es[edge_index]["bits"] = [[3,2],[2,1]]
      end
    end
    return g
  end

  function vertex_lattice(ds::Array{Int}, gates::Array, nbits; gm::Array{Int}=Array{Int}(0), st::Array{Int}=Array{Int}(0), ps::Array=Array{AbstractFloat}(0))
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
    isempty(st) || (g.vs["state"] = st)
    # g.vs["tensor"] = init_tensors(g)
    return g
  end

  function compute_energy(gate1::Function, gate2::Function, bits1::Array{Int}, bits2::Array{Int}, state1::Int, state2::Int)
    energy = 0.
    for i in 1:length(bits1)
      energy += (bin(gate1(state1), 3)[end + 1 - bits1[i]] ≠ bin(state2, 3)[end + 1 - bits2[i]])
    end
    energy
  end
  """
  Returns the energy difference for given bits, states and a new state for one vertex.
  Only one gate is needed, if the state that will be changed is the source, the gate will be applied on that
  since the state is always defined for the input.
  Use source = true if the new state is at the source of the edge, otherwise use source = false.
  """
  function compute_energydiff(gate::Function, bits1::Array{Int}, bits2::Array{Int}, state1::Int, state2::Int, state1new::Int, source::Bool)
    energy = 0.
    if source
      for i in 1:length(bits1)
        energy -= ((gate(state1)    & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state2 & 2^(bits2[i]-1)) >> (bits2[i]-1)
        energy += ((gate(state1new) & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state2 & 2^(bits2[i]-1)) >> (bits2[i]-1)
      end
    else
      for i in 1:length(bits1)
        energy -= ((gate(state1) & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state2 & 2^(bits2[i]-1)) >> (bits2[i]-1)
        energy += ((gate(state1) & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state1new & 2^(bits2[i]-1)) >> (bits2[i]-1)
      end
    end
    energy
  end

  function compute_energydifftarget(gate1::Function, gate2::Function, bits1::Array{Int}, bits2::Array{Int}, state1::Int, state2::Int, state1new::Int)
    energy = 0.
    for i in 1:length(bits1)
      energy -= (((state1)    & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state2 & 2^(bits2[i]-1)) >> (bits2[i]-1)
      energy += (((state1) & 2^(bits1[i]-1)) >> (bits1[i]-1)) $ (state1new & 2^(bits2[i]-1)) >> (bits2[i]-1)
    end
    energy
  end

  """This function tries to solve the vertex model graph as far as possible using direct computation of the gates."""
  function solve_gategraph!(g::Graph)
    statemap = g.vs["state"]
    gatemap = g.vs["gate"]
    edgelistout = edgelist(g, target=false)
    edgelistin = edgelist(g, source=false)
    edgesourcemap = edgeconnection(g, target=false)
    edgetargetmap = edgeconnection(g, source=false)
    (edgebitsinmap, edgebitsoutmap) = edgebitsmap(g, input=true, output=true)
    g.vs["state"] = solve_gategraph(g,statemap=statemap, gatemap=gatemap, edgelistin=edgelistin, edgelistout=edgelistout, edgesourcemap=edgesourcemap, edgetargetmap=edgetargetmap, edgebitsinmap=edgebitsinmap, edgebitsoutmap=edgebitsoutmap)
  end


  function solve_gategraph(g::Graph; forward::Bool=true, backward::Bool=true, statemap::Vector{Int}=Vector{Int}(), gatemap=nothing, edgelistin=nothing, edgelistout=nothing, edgesourcemap=nothing, edgetargetmap=nothing, edgebitsinmap=nothing, edgebitsoutmap=nothing)
    ngates = length(statemap)
    if forward
      for vertex in 1:ngates
        statemap[vertex] ≠ -1 && continue
        isempty(edgelistin[vertex]) && continue
        -1 in [statemap[x] for x in edgesourcemap[edgelistin[vertex]]] && continue
        statemap[vertex] = 0
        for edge in edgelistin[vertex]
          statemap[vertex] += sum(gatemap[edgesourcemap[edge]](statemap[edgesourcemap[edge]]) & 2 .^ (edgebitsoutmap[edge] .- 1) .<< (edgebitsinmap[edge] .- edgebitsoutmap[edge]))
        end
      end
    end
    if backward
      for vertex in ngates:-1:1
        statemap[vertex] ≠ -1 && continue
        isempty(edgelistout[vertex]) && continue
        -1 in [statemap[x] for x in edgetargetmap[edgelistout[vertex]]] && continue
        statemap[vertex] = 0
        for edge in edgelistout[vertex]
          statemap[vertex] += sum(statemap[edgetargetmap[edge]] & 2 .^ (edgebitsinmap[edge] .- 1) .<< (edgebitsoutmap[edge] .- edgebitsinmap[edge]))
        end
        statemap[vertex] = gatemap[vertex](statemap[vertex], true)
      end
    end
    return statemap
  end

  """ This function checks if there are mistakes in the connection between the edges of the graph."""
  function check_graph_solution(g::Graph)
    errnumber = 0
    for edge in g.es
      if g.vs[edge.source]["state"] ≡ -1 || g.vs[edge.target]["state"] ≡ -1 continue end
      if bits(g.vs[edge.target]["state"])[end+1-edge["bits"][2]] ≠ bits(g.vs[edge.source]["gate"](g.vs[edge.source]["state"]))[end+1-edge["bits"][1]]
        println("Edge $(edge.index) from gate $(edge.source) to gate $(edge.target) doesn't match.")
        errnumber += 1
      end
    end
    errnumber ≠ 0 && println("$errnumber edges don't match in the graph.")
    return
  end

  """ This function initialized the vertex model for the annealing protocol."""
  function init_annealing!(g::Graph, ds::Array{Int}, fixed_left::Array{Int}=Array{Int,2}(), fixed_right::Array{Int}=Array{Int,2}())
    g.vs[1 : end]["state"] = rand(0:7, length(g.vs["state"]))
    g.vs["fixed"] = false
    g.vs[fixed_left[:,1]]["fixed"] = true
    g.vs[end - ds[2] + fixed_right[:,1]]["fixed"] = true
    g.vs[fixed_left[:,1]]["state"] = fixed_left[:,2]
    g.vs[end - ds[2] + fixed_right[:,1]]["state"] = fixed_right[:,2]
    return
  end

  """
  Annealing function
    do_majority = true if state_counts should be returned and majority learning should be applied
    fixstates = true if the majority algorithm should fix the states
    calc_energy = true if the energy value of the system should be returned
    totalsteps is the number of steps the system should learn
    majority_steps is the number of replicas over which the system is averaged to get the majority states
    majority_cutoff is the % of the state replicas needed to fix a state
    mcsteps is the number of montecarlo steps for every annealing
    """
    function anneal!(g::Graph, ds::Array{Int};
     mcsteps::Int=2^10, majority_steps::Int=500,
     majority_cutoff::AbstractFloat=0.75, totalsteps::Int=1, calc_energy::Bool=false, fixstates::Bool=true, do_majority::Bool=true)

    do_majority && (States = SharedArray{Int}(length(g.vs), majority_steps))
    do_majority && (state_counts = Array{Array,1}(totalsteps))
    calc_energy && (energy= SharedArray{Float64}(majority_steps, totalsteps))
    for step in 1:totalsteps
      @sync @parallel for j in 1:majority_steps
        for vertex in g.vs
          vertex["fixed"] || (vertex["state"] = rand(0:7))
        end
        statemap = montecarlo!(g, mcsteps)
        do_majority && (States[:,j] = statemap)
        calc_energy && (energy[j,step] = totalunfits(g, statemap))
      end
      do_majority && (majority_steps > 1 && (state_counts[step] = majority(g, States, majority_cutoff, majority_steps, fixstates)))
    end
    calc_energy && do_majority && (return state_counts, energy)
    calc_energy && (return energy)
    do_majority && (return state_counts)
  end

  """Calculates one complete annealing from Temperature T = 1 to (almost) 0 in mcsteps montecarlo steps"""
  function montecarlo!(g::Graph, mcsteps::Int)
    gmap = g.vs["gate"]; bitsmap = g.es["bits"]
    statemap = g.vs["state"]; fixedmap = g.vs["fixed"]
    ngates = length(g.vs)
    for t in 1:mcsteps
      β = 1/(1-t/(mcsteps+1))
      for i in 1:ngates
        vertex = rand(1:ngates)
        fixedmap[vertex] && continue
        newstate = rand(0:7)
        newstate == statemap[vertex] && continue
        neighbor_vertices_index = neighbors(g, vertex)
        ΔE = 0.
        for neighbor_vertex in neighbor_vertices_index
          edge = get_eid(g, vertex, neighbor_vertex)[1]
          if g.es[edge].source == vertex
            ΔE += compute_energydiff(gmap[vertex], bitsmap[edge][1], bitsmap[edge][2], statemap[vertex], statemap[neighbor_vertex], newstate, true)
          else
            ΔE += compute_energydiff(gmap[neighbor_vertex], bitsmap[edge][1], bitsmap[edge][2], statemap[neighbor_vertex], statemap[vertex], newstate, false)
          end
        end
        (rand() < exp(- β * ΔE)) && (statemap[vertex] = newstate)
      end
    end
    statemap
  end

  function majority(g::Graph, States, majority_cutoff::AbstractFloat, majority_steps::Int, fixstates::Bool)
    state_counts = Array{Array{Float64},1}()
    for i in 1:length(g.vs)
      append!(state_counts, [counts(States[i,:], 0:7) ./ majority_steps])
      for (j, state_count) in enumerate(state_counts[i])
        if state_count ≥ majority_cutoff
          g.vs[i]["state"] = j-1
          fixstates && (g.vs[i]["fixed"] = true)
        end
      end
    end
    return state_counts
  end

  function totalunfits(g::Graph, statemap::Array{Int}=Array{Int}(0))
    isempty(statemap) && (statemap = g.vs["state"])
    energy = 0.
    for vertex in g.vs
      neighbor_vertices_index = neighbors(g, vertex, source=false, target=true)
      for neighbor_vertex_index in neighbor_vertices_index
        neighbor_vertex = g.vs[neighbor_vertex_index]
        edge = g.es[get_eid(g, vertex.index, neighbor_vertex.index)[1]]
        energy += compute_energy(vertex["gate"], neighbor_vertex["gate"], edge["bits"][1], edge["bits"][2], statemap[vertex.index], statemap[neighbor_vertex_index])
      end
    end
    energy
  end

  """
  BETA ALERT: FUNCTION SHOULD WORK PROPERLY, BUT WAS NOT YET THROUGLY TESTED.
  IF ANY ERROR SHOULD OCCUR PLESE REPORT!
  Function checks if the graph with the given parameters has a
  unique solution.
  The function returns the states of the fixed gates on each side.
  If No unique solution is found, the function tries over up to
  100 different configurations of states for of the fixed gates.
  If gm (gatemap) is given, it will try if that particular configuration
  has a unique solution.

  ds are the dimensions of the graph (i.e. (lx,ly))
  fixed_left are the states to be fixed on the left side.
  fixed_right are the states to be fixed on the right side.
  g is the optional graph that can be passed to work on
  gm is the optional gatemap
  gates is a list with the functions of the gates of the graph
  nbits are the number of bits of each gate
  """
  function unique_solution(ds::Array{Int}, fixed_left::Array{Int}, fixed_right::Array{Int}; g=nothing, gm=Array{Int}(0), gates=nothing,ps=[1,1,1,1,1], nbits::Int=3, solution::Bool=false, trials::Int=100)
    # Number of possible states per gate
    nstates = 2^nbits
    # If no graph provided, create a new one
    g ≡ nothing ? (g = vertex_lattice(ds, gates, nbits, gm=gm, ps=ps)) : gm = 1
    # gm = 1 for the return statement
    statemap = g.vs["state"]
    gatemap = g.vs["gate"]
    edgelistout = edgelist(g, target=false)
    edgelistin = edgelist(g, source=false)
    edgesourcemap = edgeconnection(g, target=false)
    edgetargetmap = edgeconnection(g, source=false)
    (edgebitsinmap, edgebitsoutmap) = edgebitsmap(g, input=true, output=true)
    statemap[:] = -1
    # Chooses which way of checking is faster (from right to left or left to right)
    if length(fixed_left) ≥ length(fixed_right)
      statemap[fixed_left] = rand(0:nstates-1, length(fixed_left))
      not_fixed_left = Array{Int,1}()
      for i in 1:ds[2]
        statemap[i] == -1 && push!(not_fixed_left, i)
      end
      unfixed = length(not_fixed_left)

      # left_states and right_states are 2d arrays with dimensions [2,ds[1]]
      # where the dimension [0,:] is the gate row_index and [1,:] the state
      left_states = hcat(1:ds[2], zeros(Int, ds[2]))
      left_states[fixed_left, 2] = statemap[fixed_left]
      statemap[left_states[:,1]] = left_states[:,2]
      initstate = rand(1:nstates^unfixed)
      pos = idxtopos(initstate, ones(Int, unfixed) * nstates) .- 1
      left_states[not_fixed_left, 2] = pos
      statemap[:] = -1
      statemap[left_states[:,1]] = left_states[:,2]
      sol = solve_gategraph(g,statemap=copy(statemap), gatemap=gatemap, edgelistin=edgelistin, edgelistout=edgelistout, edgesourcemap=edgesourcemap, edgetargetmap=edgetargetmap, edgebitsinmap=edgebitsinmap, edgebitsoutmap=edgebitsoutmap)    # Those right_states will be checked for uniqueness
      right_states = hcat(1:ds[2], sol[end - ds[2] + 1:end])

      for i in 1:nstates^unfixed
        # pos iterates over the states of the unfixed gates i.e. [0,0,1] ..
        # [0,0,2]
        i == initstate && continue
        pos = idxtopos(i, ones(Int, unfixed) * nstates) .- 1
        left_states[not_fixed_left, 2] = pos
        statemap[:] = -1
        statemap[left_states[:,1]] = left_states[:,2]
        solve_gategraph(g,statemap=statemap, gatemap=gatemap, edgelistin=edgelistin, edgelistout=edgelistout, edgesourcemap=edgesourcemap, edgetargetmap=edgetargetmap, edgebitsinmap=edgebitsinmap, edgebitsoutmap=edgebitsoutmap)
        # saved right_states before the loop are compared to the current
        # ones
        if any((right_states[fixed_right, 2] - statemap[end - ds[2] + fixed_right]) .≠ 0)
          continue
        else
        println("Found double solution")
          if trials > 0
            println("trying new random fixed states")
            if gm == Array{Int}(0)
              return unique_solution(ds, fixed_left, fixed_right, gates=gates, trials=trials-1, solution=solution, ps=ps)
            else
              return unique_solution(ds, fixed_left, fixed_right, g=g, trials=trials-1, solution=solution, ps=ps)
            end
          else
            println("still no solution")
            return
          end
        end
      end
    else
      statemap[end - ds[2] + fixed_right] = rand(0:nstates-1, length(fixed_right))
      not_fixed_right = Array{Int,1}()
      for i in 1:ds[2]
        (statemap[end-ds[2]+i] == -1) && push!(not_fixed_right, i)
      end
      unfixed = length(not_fixed_right)

      # left_states and right_states are 2d arrays with dimensions [2,ds[1]]
      # where the dimension [0,:] is the gate row_index and [1,:] the state
      right_states = hcat(1:ds[2], zeros(Int, ds[2]))
      right_states[fixed_right, 2] = statemap[end - ds[2]  + fixed_right]
      statemap[end - ds[2] + right_states[:,1]] = right_states[:,2]
      initstate = rand(1:nstates^unfixed)
      pos = idxtopos(initstate, ones(Int, unfixed) * nstates) .- 1
      right_states[not_fixed_right, 2] = pos
      statemap[:] = -1
      statemap[end - ds[2] + right_states[:,1]] = right_states[:,2]
      sol = solve_gategraph(g,statemap=copy(statemap), gatemap=gatemap, edgelistin=edgelistin, edgelistout=edgelistout, edgesourcemap=edgesourcemap, edgetargetmap=edgetargetmap, edgebitsinmap=edgebitsinmap, edgebitsoutmap=edgebitsoutmap)
      # These left_states will be checked for uniqueness
      left_states = hcat(1:ds[2], sol[1:ds[2]])

      for i in 1:nstates^unfixed
        # pos iterates over the states of the unfixed gates i.e. [0,0,1] .. [0,0,2]
        i == initstate && continue
        pos = idxtopos(i, ones(Int, unfixed) * nstates) .- 1
        right_states[not_fixed_right, 2] = pos
        statemap[:] = -1
        statemap[end - ds[2] + right_states[:,1]] = right_states[:,2]
        solve_gategraph(g,statemap=statemap, gatemap=gatemap, edgelistin=edgelistin, edgelistout=edgelistout, edgesourcemap=edgesourcemap, edgetargetmap=edgetargetmap, edgebitsinmap=edgebitsinmap, edgebitsoutmap=edgebitsoutmap)
        # saved left_states before the loop are compared to the current ones
        if (any((left_states[fixed_left, 2] - statemap[fixed_left]) .≠ 0))
          continue
        else
        println("Found double solution")
          if(trials >0)
            println("trying new random fixed states")
            if gm == Array{Int}(0)
              return unique_solution(ds, fixed_left, fixed_right, gates=gates, trials=trials-1, solution=solution, ps=ps)
            else
              return unique_solution(ds, fixed_left, fixed_right, g=g, trials=trials-1, solution=solution, ps=ps)
            end
          else
            println("still no solution")
            return
          end
        end
      end
    end
    # If no graph, nor gatemap was provided, then the gatemap
    # produced in the function is passed.
    # Returned are only the states of the fixed gates.
    if gm == Array{Int}(0)
      if solution
        return (left_states[fixed_left, 2], right_states[fixed_right, 2], g.vs["type"][:], sol)
      else
        return (left_states[fixed_left, 2], right_states[fixed_right, 2], g.vs["type"][:])
      end
    else
      if solution
        return (left_states[fixed_left, 2], right_states[fixed_right, 2], sol)
      else
        return (left_states[fixed_left, 2], right_states[fixed_right, 2])
      end
    end
  end

  function energymeasure(fname, startdecade::Int, enddecade::Int, bins::Int, ds, gates, ps)
    energy = SharedArray{Float64}(bins, (enddecade - startdecade) + 1)
    for decade in startdecade:enddecade
      @sync @parallel for j in 1:bins
        g = vertex_lattice(ds,gates, 3, ps=ps)
        g.vs["fixed"] = false
        g.vs["state"] = rand(0:7,length(g.vs))
        statemap = montecarlo!(g, 2^decade)
        energy[j,decade - startdecade + 1] = totalunfits(g, statemap)
      end
    end
    open(fname,"w") do f
      writedlm(f,energy)
    end
    return energy
  end

  function statecounts(fname, decades::Int, ds, gates, ps,totalsteps,majority_steps, fixed_left, fixed_right)
    left_states = hcat(fixed_left, zeros(Int, length(fixed_left)))
    right_states = hcat(fixed_right, zeros(Int, length(fixed_right)))
    left_states[:,2], right_states[:,2], typ, solution =unique_solution(ds, fixed_left, fixed_right,  solution=true,gates = gates,ps=ps, trials = 100)
    g = vertex_lattice(ds,gates, 3, gm=typ)
    init_annealing!(g, ds, left_states, right_states)
    state_counts, energy = anneal!(g, ds; mcsteps=2^decades, majority_steps=majority_steps, majority_cutoff=0.75,
    totalsteps=totalsteps, calc_energy=true, fixstates=true, do_majority=true)
    energy /= prod(ds)
    open(fname * ".txt","w") do f
      writedlm(f, energy)
    end
    save(fname * ".jld", "state_counts", state_counts, "gatemap", g.vs["type"], "solution", solution, "ds", ds, "fixed", g.vs["fixed"])
  end

  function get_uniquesolutions(ds, fixed_left, fixed_right; gates=[idid, toff, swsw, swid, idsw], ps=[1,1,1,1,1], num_graphs=10, fname=nothing, append_fname=false, trials=1000)
    left_states = hcat(fixed_left, zeros(Int, length(fixed_left)))
    right_states = hcat(fixed_right, zeros(Int, length(fixed_right)))
    (fname == nothing) && (fname = "unique_solution_graph_$(ds[1])x$(ds[2])_fl$(length(fixed_left))_fr$(length(fixed_right))")
    (typeof(fname) == Int) && (fname = "$fname" * "unique_solution_graph_$(ds[1])x$(ds[2])_fl$(length(fixed_left))_fr$(length(fixed_right))")
    append_fname && (fname = fname * "unique_solution_graph_$(ds[1])x$(ds[2])_fl$(length(fixed_left))_fr$(length(fixed_right))")
    for i in 1:num_graphs
      left_states[:,2], right_states[:,2], gatetype, solution = unique_solution(ds, fixed_left, fixed_right, gates=gates, ps=ps, solution=true, trials=trials)
      save(fname*"_$i.jld", "left_states", left_states, "right_states", right_states, "gatetype", gatetype, "solution", solution)
    end
  end
  end
###################################################################################################################
