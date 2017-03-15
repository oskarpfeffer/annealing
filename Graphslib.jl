
"""
Fast Graphs library for Julia with simple functionalities.
  First instance February 2017 by Oskar Pfeffer


Arguments for Edge:
index, source, target, attributes

Arguments for Vertex:
index, Set(Edges), attributes

Arguments for Graph:
  Array{Vertex}, Array{Edge}, attributes

  functions:
  expansion of setindex! and getindex for type Edge, Vertex and Graph
  newgraph
  neighbors
  add_vertices!
  add_edges!
  get_eid
  remove_edges!
  remove_vertices!
"""
module Graphslib

  export Graph, Edge, Vertex, newgraph, neighbors, add_vertices!, add_edges!, get_eid, remove_edges!, remove_vertices!
  import Base.setindex!
  import Base.getindex

  abstract GraphTypes #New supertype definition

  type Edge <: GraphTypes
    index::Int
    source::Int
    target::Int
    attributes::Dict{Any,Any}
    function Edge(index, source, target, attributes=Dict())
      new(index, source, target, attributes)
    end
  end

  type Vertex <: GraphTypes
    index::Int
    edges::IntSet
    attributes::Dict
    function Vertex(index::Int, attributes::Dict=Dict())
      new(index, IntSet(), attributes)
    end

  end

  type Graph <: GraphTypes
    vs::Array{Vertex,1}
    es::Array{Edge,1}
    attributes::Dict
    function Graph(vs::Array{Vertex,1}=Array(Vertex,0), es::Array{Edge,1}=Array(Edge,0), attributes::Dict=Dict())
      new(vs, es, attributes)
    end
  end


  function setindex!(G::GraphTypes, val::Any, key::Any)
    G.attributes[key] = val
  end

  function setindex!(G::Vector{Vertex}, val, key::String)
    for vertex in G
      vertex[key]=val
    end
  end

  function setindex!(G::Vector{Vertex}, val::Array, key::String)
    length(G) ≠ length(val) && error("length(G) ≠ length(value)")
    i = 0
    for vertex in G
      i += 1
      vertex[key] = val[i]
    end
    return G
  end

  function setindex!(G::Vector{Edge}, val, key::String)
    for edge in G
      edge[key] = val
    end
  end

  function setindex!(G::Vector{Edge}, val::Array, key::String)
    length(G) ≠ length(val) && error("length(G) ≠ length(value)")
    i = 0
    for edge in G
      i += 1
      edge[key] = val[i]
    end
  end
  function getindex(G::GraphTypes, key::String)
    return G.attributes[key]
  end
  function getindex(G::Vector{Vertex}, key::String)
    A = []
    for vertex in G
      push!(A, vertex[key])
    end
    try A = Array{typeof(A[1]),1}(A) end
    return A
  end
  function getindex(G::Vector{Edge}, key::String)
    A = []
    for edge in G
      push!(A, edge[key])
    end
    try A = Array{typeof(A[1]),1}(A) end
    return A
  end

  """ Returns a Graph with N vertices."""
  function newgraph(NVertex::Int=0)
    Vertices = Array(Vertex,NVertex)
    for i in 1:NVertex
      Vertices[i] = Vertex(i)
    end
    Graph(Vertices)
  end

  """
  Returns the neighbors of a given vertex in a given graph.
  If source is true it returns the neighbors connected through a graph where
  they are the source.
  If target is true it returns the neighbors connected
  through a graph where they are the target.
  If both are true, it returns all the neighbors.
  """
  function neighbors(g::Graph, vertex::Vertex, source::Bool=true, target::Bool=true)
    neighbors = Array(Int,0)
    for edge in vertex.edges
      if g.es[edge].source == vertex.index
        if target
          append!(neighbors, g.es[edge].target)
        end
      else
        if source
          append!(neighbors, g.es[edge].source)
        end
      end
    end
    return neighbors
  end

  function neighbors(g::Graph, vertex::Int, source::Bool=true, target::Bool=true)
    return neighbors(g, g.vs[vertex], source, target)
  end

  """
  Adds N vertices to the graph g with optional given attributes.
  Attributes are the same for every vertex.
  """
  function add_vertices!(g::Graph, N::Int, attributes::Dict=Dict())
    for i in 1:N
      push!(g.vs, Vertex(length(g.vs)+1, attributes))
    end
  end

  """
  Adds an edge from vertex source to vertex target with given attributes in graph g.
  """
  function add_edges!(g::Graph, source, target, attributes::Dict=Dict())
    source==target && error("Cannot connect a vertex to himself.")
    source > length(g.vs) && error("Given source vertex does not exist.")
    target > length(g.vs) && error("Given target vertex does not exist.")
    push!(g.es, Edge(length(g.es)+1, source, target, attributes))
    push!(g.vs[source].edges, length(g.es))
    push!(g.vs[target].edges,length(g.es))
    return nothing
  end

  """Adds multiple edges with same attribute for every edge."""
  function add_edges!(g::Graph, vertices::Array{Int, 2}, attributes::Dict=Dict())
    for edge in vertices
      add_edges!(g, edge..., attributes)
    end
  end

  """Adds multiple edges with different attribute(Array{Dict}) for every edge."""
  function add_edges!(g::Graph, vertices::Array{Int, 2}, attributes::Array{Dict})
    for i in 1:length(vertices)
      add_edges!(g, vertices[i,:]..., attributes[i])
    end
  end

  """
  Returns an Array{Int,1} containing all the edges connecting
  two vertices in the graph g.
  """
  function get_eids(g::Graph, source::Int, target::Int)
    indices = Array{Int,1}()
    for edge in g.es
      if ((edge.source == source && edge.target == target) || (edge.source == target && edge.target == source) )
        push!(indices, edge.index)
      end
    end
    return indices
  end
  function get_eid(g::Graph, source::Int, target::Int)
    for edge in g.es
      if ((edge.source == source && edge.target == target) || (edge.source == target && edge.target == source) )
        return edge.index
      end
    end
  end
  """
  Returns an Array{Int,1} containing all the edges connecting
  two vertices in the graph g.
  """
  function get_eids(g::Graph, source::Vertex, target::Vertex)
    return get_eids(g, source.index, target.index)
  end
  function get_eid(g::Graph, source::Vertex, target::Vertex)
    return get_eid(g, source.index, target.index)
  end

  """Removes edge in graph g by its index."""
  function remove_edges!(g::Graph, index::Int)
    pop!(g.vs[g.es[index].source].edges, index)
    pop!(g.vs[g.es[index].target].edges, index)
    for edge in g.es[index+1:end]
      pop!(g.vs[g.es[edge.index].source].edges, edge.index)
      pop!(g.vs[g.es[edge.index].target].edges, edge.index)
      push!(g.vs[g.es[edge.index].source].edges, edge.index-1)
      push!(g.vs[g.es[edge.index].target].edges, edge.index-1)
      edge.index -= 1
    end
    deleteat!(g.es, index)
    return nothing
  end

  """Removes multiple edges in graph g by their index(Vector{Int})."""
  function remove_edges!(g::Graph, index::Vector{Int})
    index = sort(unique(index))
    for i in index[end:-1:1]
      remove_edges!(g, i)
    end
  end

  """Removes multiple edges in graph g by their index(Vector{Int})."""
  function remove_edges!(g::Graph, index::IntSet)
    index_copy = copy(index)
    while length(index_copy)>0
      remove_edges!(g, maximum(index_copy))
      pop!(index_copy, maximum(index_copy))
    end
  end

  """Removes given edge in graph g."""
  function remove_edges!(g::Graph, edge::Edge)
    remove_edges!(g, edge.index)
  end

  """Removes given edges in graph g."""
  function remove_edges!(g::Graph, edges::Vector{Edge})
    to_remove = Vector{Int}()
    for edge in edges
      push!(to_remove, edge.index)
    end
    remove_edges!(g, to_remove)
  end

  """Removes vertex in graph g by its index."""
  function remove_vertices!(g::Graph, index::Int)
    remove_edges!(g, g.vs[index].edges)
    for vertex in g.vs[index+1:end]
      for edge_index in vertex.edges
        edge = g.es[edge_index]
        if edge.source == vertex.index
          edge.source-=1
        elseif edge.target == vertex.index
          edge.target -= 1
        end
      end
      vertex.index-=1
    end
    deleteat!(g.vs, index)
  end

  """Removes vertices in graph g by their indices."""
  function remove_vertices!(g::Graph, index::Vector{Int})
    index = Set(index)
    remove_vertices!(g, index)
  end

  """Removes vertices in graph g by their indices."""
  function remove_vertices!(g::Graph, index::IntSet)
    while length(index)>0
      remove_vertices!(g, maximum(index))
      pop!(index, maximum(index))
    end
  end

  """Removes given vertex in graph g."""
  function remove_vertices!(g::Graph, vertex::Vertex)
    remove_vertices!(g, vertex.index)
  end

  """Removes given vertices in graph g."""
  function remove_vertices!(g::Graph, vertices::Vector{Vertex})
    to_remove = Vector{Int}()
    for vertex in vertices
      push!(to_remove, vertex.index)
    end
    remove_vertices!(g, to_remove)
  end

end
