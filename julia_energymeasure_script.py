ps = ["[0.25  , 0.    , 0.25  , 0.25  , 0.25  ]",
      "[0.2375, 0.05  , 0.2375, 0.2375, 0.2375]",
      "[0.225 , 0.1   , 0.225 , 0.225 , 0.225 ]",
      "[0.2125, 0.15  , 0.2125, 0.2125, 0.2125]",
      "[0.2   , 0.2   , 0.2   , 0.2   , 0.2   ]",
      "[0.175 , 0.3   , 0.175 , 0.175 , 0.175 ]",
      "[0.15  , 0.4   , 0.15  , 0.15  , 0.15  ]",
      "[0.125 , 0.5   , 0.125 , 0.125 , 0.125 ]",
      "[0.1   , 0.6   , 0.1   , 0.1   , 0.1   ]",
      "[0.05  , 0.8   , 0.05  , 0.05  , 0.05  ]",
      "[0.    , 1.    , 0.    , 0.    , 0.    ]"]
for (j,PS) in enumerate(ps):
    percent = [0,5,10,15,20,30,40,50,60,80,100]
    for iteration in [1]:
        for L in [5, 10, 20]:
            if L ==5:
                
                for (i,W) in enumerate([50, 75, 100, 200, 300]):
                    file = open("startenergy{}_{}x{}_{}.jl".format(percent[j],L,W, iteration), "w")
                    file.write("push!(LOAD_PATH, pwd())\n")
                    file.write("using annealing\n")
                    file.write("gates = [idid, toff, swsw, swid, idsw]\n")
                    replicas = [200, 150, 100, 50, 32]
                    file.write("ds = [{}, {}]\n".format(W,L))
                    file.write("ps = "+PS+"\n")
                    file.write("energymeasure(\"compile.txt\", 1, 1, 1, ds, gates, ps)\n")
                    file.write("@time energymeasure(\"dat/{}x{}_{}toff_20dec_{}replicas_{}.txt\", 10, 20, {}, ds, gates, ps)".format(W, L, percent[j], replicas[i], iteration, replicas[i]))
            if L == 10:
                for (i,W) in enumerate([50, 100, 200]):
                    file = open("startenergy{}_{}x{}_{}.jl".format(percent[j],L,W, iteration), "w")
                    file.write("push!(LOAD_PATH, pwd())\n")
                    file.write("using annealing\n")
                    file.write("gates = [idid, toff, swsw, swid, idsw]\n")
                    replicas = [100, 50, 32]
                    file.write("ds = [{}, {}]\n".format(W,L))
                    file.write("ps = "+PS+"\n")
                    file.write("energymeasure(\"compile.txt\", 1, 1, 1, ds, gates, ps)\n")
                    file.write("@time energymeasure(\"dat/{}x{}_{}toff_20dec_{}replicas_{}.txt\", 10, 20, {}, ds, gates, ps)".format(W,L,percent[j],replicas[i],iteration, replicas[i]))
            if L == 20:
                for (i,W) in enumerate([50, 100, 200]):
                    file = open("startenergy{}_{}x{}_{}.jl".format(percent[j],L,W, iteration), "w")
                    file.write("push!(LOAD_PATH, pwd())\n")
                    file.write("using annealing\n")
                    file.write("gates = [idid, toff, swsw, swid, idsw]\n")
                    replicas = [64, 32, 16]
                    file.write("ds = [{}, {}]\n".format(W,L))
                    file.write("ps = "+PS+"\n")
                    file.write("energymeasure(\"compile.txt\", 1, 1, 1, ds, gates, ps)\n")
                    file.write("@time energymeasure(\"dat/{}x{}_{}toff_20dec_{}replicas_{}.txt\", 10, 20, {}, ds, gates, ps)".format(W,L,percent[j],replicas[i],iteration, replicas[i]))
