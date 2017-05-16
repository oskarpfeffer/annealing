
module toymodel
  export init_toymodel, propagation, fixbits, convert_to_state
  using Graphslib



  function init_toymodel(g, ds,fixed_left=[], fixed_right=[])
    for vertex in g.vs
      vertex["bitprobability"] = [[0.5, 0.5, 0.5], [0.5, 0.5, 0.5]]
      vertex["fixed"] = [0,0,0]
    end
    gate_map = g.vs["type"]
    bitprob = g.vs["bitprobability"]
    for i in fixed_left
      st = g.vs[i]["state"]
      bitprob[i][1] = [(st&1), (st&2)>>1, (st&4)>>2]
      g.vs[i]["fixed"] = [1,1,1]
      if gate_map[i] == 1
        bitprob[i][2] = bitprob[i][1]
      elseif gate_map[i] == 2
        bitprob[i][2][1] = bitprob[i][1][1]
        bitprob[i][2][2] = bitprob[i][1][2]
        bitprob[i][2][3] = bitprob[i][1][3] * (1 - bitprob[i][1][1] * bitprob[i][1][2]) + (1 - bitprob[i][1][3]) * bitprob[i][1][1] * bitprob[i][1][2]
      elseif gate_map[i] == 3
        bitprob[i][2][1] = bitprob[i][1][3]
        bitprob[i][2][2] = bitprob[i][1][1]
        bitprob[i][2][3] = bitprob[i][1][2]
      elseif gate_map[i] == 4
        bitprob[i][2][1] = bitprob[i][1][1]
        bitprob[i][2][2] = bitprob[i][1][3]
        bitprob[i][2][3] = bitprob[i][1][2]
      elseif gate_map[i] == 5
        bitprob[i][2][1] = bitprob[i][1][2]
        bitprob[i][2][2] = bitprob[i][1][1]
        bitprob[i][2][3] = bitprob[i][1][3]
      end
    end
    for i in fixed_right
      index = length(g.vs) - ds[2] + i
      st = g.vs[index]["state"]
      bitprob[index][1] = [(st&1), (st&2)>>1, (st&4)>>2]
      g.vs[index]["fixed"] = [1,1,1]
      if gate_map[index] == 1
        bitprob[index][2] = bitprob[index][1]
      elseif gate_map[index] == 2
        bitprob[index][2][1] = bitprob[index][1][1]
        bitprob[index][2][2] = bitprob[index][1][2]
        bitprob[index][2][3] = bitprob[index][1][3] * (1 - bitprob[index][1][1] * bitprob[index][1][2]) + (1 - bitprob[index][1][3]) * bitprob[index][1][1] * bitprob[index][1][2]
      elseif gate_map[index] == 3
        bitprob[index][2][1] = bitprob[index][1][3]
        bitprob[index][2][2] = bitprob[index][1][1]
        bitprob[index][2][3] = bitprob[index][1][2]
      elseif gate_map[index] == 4
        bitprob[index][2][1] = bitprob[index][1][1]
        bitprob[index][2][2] = bitprob[index][1][3]
        bitprob[index][2][3] = bitprob[index][1][2]
      elseif gate_map[index] == 5
        bitprob[index][2][1] = bitprob[index][1][2]
        bitprob[index][2][2] = bitprob[index][1][1]
        bitprob[index][2][3] = bitprob[index][1][3]
      end
    end
  end

  function propagation(g, loss, right=true,left=true)
    loss = 1 - loss
    bitprob = g.vs["bitprobability"]
    fixed_map = g.vs["fixed"]
    gate_map = g.vs["type"]
    if right
      for vertex in 1:length(g.vs)
        sources_vertices = neighbors(g, vertex, source=true, target=false)
        for source in sources_vertices
          edge = g.es[get_eid(g, vertex, source)]
          for i in 1:length(edge["bits"][1])
            fixed_map[vertex][edge["bits"][2][i]] == 1 && continue
            bitprob[vertex][1][edge["bits"][2][i]] = (loss * bitprob[source][2][edge["bits"][1][i]]) + (0.5 * (1-loss))
            if gate_map[vertex] == 1
            end
          end
        end
        if gate_map[vertex] == 1
          bitprob[vertex][2][1] = bitprob[vertex][1][1]
          bitprob[vertex][2][2] = bitprob[vertex][1][2]
          bitprob[vertex][2][3] = bitprob[vertex][1][3]
        elseif gate_map[vertex] == 2
          bitprob[vertex][2][1] = bitprob[vertex][1][1]
          bitprob[vertex][2][2] = bitprob[vertex][1][2]
          bitprob[vertex][2][3] = bitprob[vertex][1][3] * (1 - bitprob[vertex][1][1] * bitprob[vertex][1][2]) + (1 - bitprob[vertex][1][3]) * bitprob[vertex][1][1] * bitprob[vertex][1][2]
        elseif gate_map[vertex] == 3
          bitprob[vertex][2][1] = bitprob[vertex][1][3]
          bitprob[vertex][2][2] = bitprob[vertex][1][1]
          bitprob[vertex][2][3] = bitprob[vertex][1][2]
        elseif gate_map[vertex] == 4
          bitprob[vertex][2][1] = bitprob[vertex][1][1]
          bitprob[vertex][2][2] = bitprob[vertex][1][3]
          bitprob[vertex][2][3] = bitprob[vertex][1][2]
        elseif gate_map[vertex] == 5
          bitprob[vertex][2][1] = bitprob[vertex][1][2]
          bitprob[vertex][2][2] = bitprob[vertex][1][1]
          bitprob[vertex][2][3] = bitprob[vertex][1][3]

        # if gate_map[vertex] == 2
        #   bitprob[vertex][1][1] = bitprob[vertex][2][1] = (bitprob[vertex][1][1] + bitprob[vertex][2][1]) * 0.5
        #   bitprob[vertex][1][2] = bitprob[vertex][2][2] = (bitprob[vertex][1][2] + bitprob[vertex][2][2]) * 0.5
        #   bitprob[vertex][1][3] = bitprob[vertex][2][3] = (bitprob[vertex][1][3] + bitprob[vertex][2][3]) * 0.5
        # elseif gate_map[vertex] == 1
        #   bitprob[vertex][1][1] = bitprob[vertex][2][1] = (bitprob[vertex][1][1] + bitprob[vertex][2][1]) * 0.5
        #   bitprob[vertex][1][2] = bitprob[vertex][2][2] = (bitprob[vertex][1][2] + bitprob[vertex][2][2]) * 0.5
        #   bitprob[vertex][1][3] = (bitprob[vertex][1][3] + bitprob[vertex][2][3] * (1 - bitprob[vertex][2][1] * bitprob[vertex][2][2]) + (1 - bitprob[vertex][2][3]) * bitprob[vertex][2][1] * bitprob[vertex][2][2] ) * 0.5
        #   bitprob[vertex][2][3] = bitprob[vertex][1][3] * (1 - bitprob[vertex][1][1] * bitprob[vertex][1][2]) + (1 - bitprob[vertex][1][3]) * bitprob[vertex][1][1] * bitprob[vertex][1][2]
        end
      end
    end
    if left
      for vertex in length(g.vs):-1:1
        target_vertices = neighbors(g, vertex, source=false, target=true)
        for target in target_vertices
          edge = g.es[get_eid(g, vertex, target)]
          for i in 1:length(edge["bits"][1])
            fixed_map[vertex][edge["bits"][1][i]] == 1 && continue
            bitprob[vertex][2][edge["bits"][1][i]] = (loss * bitprob[target][1][edge["bits"][2][i]]) + (0.5 * (1-loss))
          end
        end
        if gate_map[vertex] == 1
          bitprob[vertex][1][1] = bitprob[vertex][2][1]
          bitprob[vertex][1][2] = bitprob[vertex][2][2]
          bitprob[vertex][1][3] = bitprob[vertex][2][3]
        elseif gate_map[vertex] == 2
          bitprob[vertex][1][1] = bitprob[vertex][2][1]
          bitprob[vertex][1][2] = bitprob[vertex][2][2]
          bitprob[vertex][1][3] = bitprob[vertex][2][3] * (1 - bitprob[vertex][2][1] * bitprob[vertex][2][2]) + (1 - bitprob[vertex][2][3]) * bitprob[vertex][2][1] * bitprob[vertex][2][2]
        elseif gate_map[vertex] == 3
          bitprob[vertex][1][1] = bitprob[vertex][2][2]
          bitprob[vertex][1][2] = bitprob[vertex][2][3]
          bitprob[vertex][1][3] = bitprob[vertex][2][1]
        elseif gate_map[vertex] == 4
          bitprob[vertex][1][1] = bitprob[vertex][2][1]
          bitprob[vertex][1][2] = bitprob[vertex][2][3]
          bitprob[vertex][1][3] = bitprob[vertex][2][2]
        elseif gate_map[vertex] == 5
          bitprob[vertex][1][1] = bitprob[vertex][2][2]
          bitprob[vertex][1][2] = bitprob[vertex][2][1]
          bitprob[vertex][1][3] = bitprob[vertex][2][3]
        # if gate_map[vertex] == 2
        #   bitprob[vertex][1][1] = bitprob[vertex][2][1] = (bitprob[vertex][1][1] + bitprob[vertex][2][1]) * 0.5
        #   bitprob[vertex][1][2] = bitprob[vertex][2][2] = (bitprob[vertex][1][2] + bitprob[vertex][2][2]) * 0.5
        #   bitprob[vertex][1][3] = bitprob[vertex][2][3] = (bitprob[vertex][1][3] + bitprob[vertex][2][3]) * 0.5
        # elseif gate_map[vertex] == 1
        #   bitprob[vertex][1][1] = bitprob[vertex][2][1] = (bitprob[vertex][1][1] + bitprob[vertex][2][1]) * 0.5
        #   bitprob[vertex][1][2] = bitprob[vertex][2][2] = (bitprob[vertex][1][2] + bitprob[vertex][2][2]) * 0.5
        #   bitprob[vertex][1][3] = (bitprob[vertex][1][3] + bitprob[vertex][2][3] * (1 - bitprob[vertex][2][1] * bitprob[vertex][2][2]) + (1 - bitprob[vertex][2][3]) * bitprob[vertex][2][1] * bitprob[vertex][2][2] ) * 0.5
        #   bitprob[vertex][2][3] = bitprob[vertex][1][3] * (1 - bitprob[vertex][1][1] * bitprob[vertex][1][2]) + (1 - bitprob[vertex][1][3]) * bitprob[vertex][1][1] * bitprob[vertex][1][2]
        end
      end
    end
    g.vs["bitprobability"] = bitprob
  end


  function convert_to_state(g, cutoff=0.25)
    bitprob = deepcopy(g.vs["bitprobability"])

    for vertex in 1:length(g.vs)
      st = 0
      for i in 1:length(bitprob[vertex][1])
        if bitprob[vertex][1][i] < cutoff
          bitprob[vertex][1][i] = 0
        elseif bitprob[vertex][1][i] > 1-cutoff
          bitprob[vertex][1][i] = 1
        else
          bitprob[vertex][1][i] = -1
          st = -1
        end
        st == -1 && continue
        st += round(Int, bitprob[vertex][1][i]) << (i-1)
      end
      g.vs[vertex]["state"]=st
    end
  end
  function fixbits(g, cutoff=0.1)
    bitprob = g.vs["bitprobability"]
    fixed_map = g.vs["fixed"]
    gate_map = g.vs["type"]

    for vertex in 1:length(g.vs)
      for i in 1:length(bitprob[vertex][1])
        if bitprob[vertex][1][i] < cutoff
          bitprob[vertex][1][i] = 0.
          fixed_map[vertex][i] = 1
        elseif bitprob[vertex][1][i] > 1-cutoff
          bitprob[vertex][1][i] = 1.
          fixed_map[vertex][i] = 1
        end
      end
      if gate_map == 1
        bitprob[vertex][2] = bitprob[vertex][1]
      elseif gate_map == 2
        bitprob[vertex][2][1] = bitprob[vertex][1][1]
        bitprob[vertex][2][2] = bitprob[vertex][1][2]
        bitprob[vertex][2][3] = bitprob[vertex][1][3] * (1 - bitprob[vertex][1][1] * bitprob[vertex][1][2]) + (1 - bitprob[vertex][1][3]) * bitprob[vertex][1][1] * bitprob[vertex][1][2]
      elseif gate_map[vertex] == 3
        bitprob[vertex][2][1] = bitprob[vertex][1][3]
        bitprob[vertex][2][2] = bitprob[vertex][1][1]
        bitprob[vertex][2][3] = bitprob[vertex][1][2]
      elseif gate_map[vertex] == 4
        bitprob[vertex][2][1] = bitprob[vertex][1][1]
        bitprob[vertex][2][2] = bitprob[vertex][1][3]
        bitprob[vertex][2][3] = bitprob[vertex][1][2]
      elseif gate_map[vertex] == 5
        bitprob[vertex][2][1] = bitprob[vertex][1][2]
        bitprob[vertex][2][2] = bitprob[vertex][1][1]
        bitprob[vertex][2][3] = bitprob[vertex][1][3]
      end
    end
    g.vs["fixed"] = fixed_map
    g.vs["bitprobability"] = bitprob
  end
  function stabilize(g)
    bitprob = g.vs["bitprobability"]
    fixed_map = g.vs["fixed"]
    gate_map = g.vs["type"]
    for vertex in 1:length(g.vs)
      if gate_map[vertex] == 2
        bitprob[vertex][2][1] = bitprob[vertex][1][1]
        bitprob[vertex][2][2] = bitprob[vertex][1][2]
        bitprob[vertex][2][3] = bitprob[vertex][1][3]
      elseif gate_map[vertex] == 1
        bitprob[vertex][2][1] = bitprob[vertex][1][1]
        bitprob[vertex][2][2] = bitprob[vertex][1][2]
        bitprob[vertex][2][3] = bitprob[vertex][1][3] * (1 - bitprob[vertex][1][1] * bitprob[vertex][1][2]) + (1 - bitprob[vertex][1][3]) * bitprob[vertex][1][1] * bitprob[vertex][1][2]
      end
    end
  end
end
