module Graphslib
abstract GraphTypes
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
  function Vertex(index::Int, attributes::Dict)
    new(index, IntSet(), attributes)
  end
  function Vertex(index::Int)
    new(index, IntSet(), Dict())
  end
end

type Graph <: GraphTypes
  vs::Array{Vertex,1}
  es::Array{Edge,1}
  attributes::Dict
  function Graph(vs=Array(Vertex,0), es=Array(Edge,0), attributes=Dict())
    new(vs, es, attributes)
  end
end

import Base.setindex!
import Base.getindex
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
    append!(A, vertex[key])
  end
  return A
end
function getindex(G::Vector{Edge}, key::String)
  A = []
  for edge in G
    append!(A, edge[key])
  end
  return A
end

function newgraph(NVertex::Int=0)
  Vertices = Array(Vertex,NVertex)
  for i in 1:NVertex
    Vertices[i] = Vertex(i)
  end
  Graph(Vertices)
end


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

function add_vertices!(g::Graph, amount::Int, attributes::Dict=Dict())
  for i in 1:amount
    print(typeof(attributes))
    push!(g.vs, Vertex(length(g.vs)+1, attributes))
  end
end

function add_edges!(g::Graph, source, target, attributes::Dict=Dict())
  source==target && error("Cannot connect a vertex to himself")
  push!(g.es, Edge(length(g.es)+1, source, target, attributes))
  push!(g.vs[source].edges,length(g.es))
  push!(g.vs[target].edges,length(g.es))
  return g.es[end]
end

function add_edges!(g::Graph, vertices::Array{Int}, attributes::Dict=Dict())
  for edge in vertices
    add_edges!(g, edge..., attributes)
  end
end

function add_edges!(g::Graph, vertices::Array{Int}, attributes::Array{Dict})
  for i in 1:length(vertices)
    add_edges!(g, vertices[i]..., attributes[i])
  end
end

function get_eid(g::Graph, source::Int, target::Int)
  indices = []
  for edge in g.es
    if ((edge.source == source && edge.target == target) || (edge.source == target && edge.target == source) )
      push!(indices, edge.index)
    end
  end
  return indices
end

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
end

function remove_edges!(g::Graph, index::Vector{Int})
  index = sort(unique(index))
  for i in index[end:-1:1]
    remove_edges!(g, i)
  end
end

function remove_edges!(g::Graph, index::IntSet)
  index_copy = copy(index)
  print(index_copy)
  while length(index_copy)>0
    remove_edges!(g, maximum(index_copy))
    pop!(index_copy, maximum(index_copy))
  end
end


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

end
