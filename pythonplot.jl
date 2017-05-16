module pythonplot
  using PyCall
  using Graphslib
  using annealing
  using JLD
  unshift!(PyVector(pyimport("sys")["path"]), "")
  @pyimport igraph
  @pyimport graphplotting as gp
  @pyimport matplotlib.cm as cmap
  export coloring_annealing, coloring_toymodel, plotit, plotfromjld

  function coloring_annealing(g, solution, state_counts)
    cm = cmap.get_cmap("hot")
    for (i, state) in enumerate(solution)
      g.vs[i]["color"] = cm(state_counts[i][state+1])
    end
    return
  end

  function coloring_toymodel(g)
    cm = cmap.get_cmap("hot")
    for i in 1:length(g.vs)
      g.vs[i]["color"] = cm(sum(g.vs[i]["fixed"])/3.)
    end
  end


  function plotit(g_in; outputfile::String="graphplot.png", plot_states::Bool=true, plot_gatetypes=false)
    g = igraph.Graph()
    g[:add_vertices](length(g_in.vs))
    for edge in g_in.es
      g[:add_edge](edge.source-1, edge.target-1)
    end
    gp.setindex_vs(g, "state", g_in.vs["state"])
    gp.setindex_vs(g, "x", g_in.vs["x"])
    gp.setindex_vs(g, "y", g_in.vs["y"])
    gp.setindex_es(g, "bits", g_in.es["bits"])
    gp.setindex_vs(g, "gate", g_in.vs["gate"])
    gp.setindex_vs(g, "type", g_in.vs["type"])
    try gp.setindex_vs(g, "color", g_in.vs["color"])
    end
    gp.graphplotting(g, outputfile, plot_states, plot_gatetypes)
  end

  function plotfromjld(fname, gates = [idid, toff, swsw, swid, idsw],plot_gatetypes=false)
    a = load(fname*".jld")
    solution = a["solution"]
    ds = a["ds"]
    state_counts = a["state_counts"]
    g = vertex_lattice(ds,gates,3)
    for i in 1:length(state_counts)
      coloring(g, solution, state_counts[i])
      plotit(g, plot_states=false, outputfile = fname*"$i.png", plot_gatetypes=plot_gatetypes)
    end
  end
end
