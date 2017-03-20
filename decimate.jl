"""

    Function decimating the vertex_lattice. It removes the 3 bit SWAP-SWAP,
    ID-ID, ID-SWAP, SWAP-ID gates and reconnects the vertices consistently with
    the deleted gate.

    First instance January 2017 by Oskar Pfeffer

"""

module decimate
  push!(LOAD_PATH, pwd())
  using Graphslib
  using annealing
  export decimation
  function decimation(g::Graph)
    """
    This function will iterate over all the vertices and replace all idid, swsw,
    swid and idsw gates reconnecting the edges accordingly.
    g :: Graph

    Function doesn't work correctly if there are some neighbouring vertices
    connected by more than 1 edge. If that is the case call the function
    "bondsglue" first
    """

    # Loop over all vertices. Uses while because the variable
    # current_vertex_index doesn't change if a gate is removed.
    current_vertex_index = 1
    while current_vertex_index < length(g.vs)

      # Assigns the current vertex that is iterated over
      vertex = g.vs[current_vertex_index]

      # Toffoli gates won't be decimated in this function
      (vertex["gate"] == toff) && (current_vertex_index += 1; continue)

      # Variables will contain the edges to the left (inedges) and the edges
      # to the right (outedges) of the current vertex
      inedges = [] ; outedges = []
      # if neighbor_vertex_index < current_vertex_index then it's to its
      # left.
      for neighbor_vertex_index in neighbors(g, current_vertex_index)
        if neighbor_vertex_index < current_vertex_index
          push!(inedges, g.es[get_eid(g, neighbor_vertex_index, current_vertex_index)])
        elseif neighbor_vertex_index > current_vertex_index
          push!(outedges, g.es[get_eid(g, neighbor_vertex_index, current_vertex_index)])
        end
      end
      # If no neighbors to the left or to the right means that it's a boundary
      # vertex and won't be decimated
      (inedges == [] || outedges == []) && (current_vertex_index += 1; continue)

      ##################### REMOVE VERTICES #####################
      # Removes the vertex and connects correspondingly
      if (vertex["gate"] == idid)
        remove_idid(g, inedges, outedges)
        remove_vertices!(g, current_vertex_index)
      elseif (vertex["gate"] == swid)
        remove_swid(g, inedges, outedges)
        remove_vertices!(g, current_vertex_index)
      elseif (vertex["gate"] == idsw)
        remove_idsw(g, inedges, outedges)
        remove_vertices!(g, current_vertex_index)
      elseif (vertex["gate"] == swsw)
        remove_swsw(g, inedges, outedges)
        remove_vertices!(g, current_vertex_index)
      end
    end
    return nothing
  end

  # REMOVE FUNCTIONS
  # The methods that remove gates are all build up the same way, the only
  # difference is the condition for how the entering bits are reconnected to the
  # next vertex depending on the internal structure of the gate.
  # Loop over all the inbits and then you loop over all the outbits and check
  # every single inbit and outbit individually. Inbit and outbit are the bits that
  # connect to the vertex. I.e. inbit==0 is the upper left bit of the vertex,
  # outbit==1 is the right middle bit of the vertex.
  #
  # inedge_outbit (outedge_inbit) is needed to write the right argument "bits" for
  # the new edge. It is the Index of the ith element of the array
  # edge["bits"][0 (or 1)], not the actual value of the ith elemnet of the
  # argument "bits"


  function remove_idid(g::Graph, inedges, outedges)
    """
    Removes the 3-bit ID-ID gate and reconnects from the top:
        a -> a
        b -> b
        c -> c
    The information of the vertex to be replaced is not needed since it is
    already in the inedges and outedges (i.e. inedge.target is the current
    vertex)
    g :: Graph where the gates have to be decimated
    inedges :: incoming edges of the vertex
    outedges :: outgoing edges of the vertex
    """
    for inedge in inedges
      inedge_outbit = 1
      for inbit in inedge["bits"][2]
        for outedge in outedges
          outedge_inbit = 1
          for outbit in outedge["bits"][1]
            if ((inbit == 1 && outbit == 1) ||
                (inbit == 2 && outbit == 2) ||
                (inbit == 3 && outbit == 3))
              add_edge(g, inedge, outedge, inedge_outbit,
                       outedge_inbit)
            end
            outedge_inbit += 1
          end
        end
        inedge_outbit += 1
      end
    end
  end


  function remove_idsw(g::Graph, inedges, outedges)
    """
    Removes the 3-bit ID-SW gate and reconnects from the top:
        a -> b
        b -> a
        c -> c
    The information of the vertex to be replaced is not needed since it is
    already in the inedges and outedges (i.e. inedge.target is the current
    vertex)
    g :: Graph where the gates have to be decimated
    inedges :: incoming edges of the vertex
    outedges :: outgoing edges of the vertex
    """
    for inedge in inedges
      inedge_outbit = 1
      for inbit in inedge["bits"][2]
        for outedge in outedges
          outedge_inbit = 1
          for outbit in outedge["bits"][1]
            if ((inbit == 1 && outbit == 2) ||
                (inbit == 2 && outbit == 1) ||
                (inbit == 3 && outbit == 3))
              add_edge(g, inedge, outedge, inedge_outbit,
                       outedge_inbit)
            end
            outedge_inbit += 1
          end
        end
        inedge_outbit += 1
      end
    end
  end


  function remove_swid(g::Graph, inedges, outedges)
    """
    Removes the 3-bit SW-ID gate and reconnects from the top:
        a -> a
        b -> c
        c -> b
    The information of the vertex to be replaced is not needed since it is
    already in the inedges and outedges (i.e. inedge.target is the current
    vertex)
    g :: Graph where the gates have to be decimated
    inedges :: incoming edges of the vertex
    outedges :: outgoing edges of the vertex
    """
    for inedge in inedges
      inedge_outbit = 1
      for inbit in inedge["bits"][2]
        for outedge in outedges
          outedge_inbit = 1
          for outbit in outedge["bits"][1]
            if ((inbit == 1 && outbit == 1) ||
                (inbit == 2 && outbit == 3) ||
                (inbit == 3 && outbit == 2))
              add_edge(g, inedge, outedge, inedge_outbit,
                         outedge_inbit)
            end
            outedge_inbit += 1
          end
        end
        inedge_outbit += 1
      end
    end
  end

  function remove_swsw(g::Graph, inedges, outedges)
    """
    Removes the 3-bit SW-SW gate and reconnects from the top:
        a -> b
        b -> c
        c -> a
    The information of the vertex to be replaced is not needed since it is
    already in the inedges and outedges (i.e. inedge.target is the current
    vertex)
    g :: Graph where the gates have to be decimated
    inedges :: incoming edges of the vertex
    outedges :: outgoing edges of the vertex
    """
    for inedge in inedges
      inedge_outbit = 1
      for inbit in inedge["bits"][2]
        for outedge in outedges
          outedge_inbit = 1
          for outbit in outedge["bits"][1]
            if ((inbit == 1 && outbit == 2) ||
                (inbit == 2 && outbit == 3) ||
                (inbit == 3 && outbit == 1))
              add_edge(g, inedge, outedge, inedge_outbit,
                         outedge_inbit)
            end
            outedge_inbit += 1
          end
        end
        inedge_outbit += 1
      end
    end
  end

  function add_edge(g::Graph, inedge, outedge, inedge_outbit::Int, outedge_inbit::Int)
    """
    Function adds a new edge connecting inedge.source to outedge.target and
    if there is already an edge, it only adds the "bits" argument to the
    existing edge.
    g :: Graph in which edge has to be add_edge
    inedge :: input_edge of the removed vertex
    outedge :: output_edge of the removed vertex
    inedge_outbit :: Index of the bit in the inedge that will be reconnected
    outedge_outbit :: Index of the bit in the outedge that will be reconnected
    """
    # Check if there is already an edge connecting the vertices and returns
    # the edge index or -1 if there is none
    edge0 = get_eid(g, inedge.source, outedge.target)

    # No edge already existing
    if (edge0 == nothing)
      add_edges!(g, inedge.source, outedge.target)
      g.es[length(g.es)]["bits"] =  [[inedge["bits"][1][inedge_outbit]],
                                    [outedge["bits"][2][outedge_inbit]]]
    # Edge exists already, append the new "bits" to the existing edge
    else
      edge0 = g.es[edge0]
      edge0["bits"] = [append!(edge0["bits"][1], inedge["bits"][1][inedge_outbit]),
           append!(edge0["bits"][2], outedge["bits"][2][outedge_inbit])]
    end
  end
end
