
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
                    file = open("job_startenergy{}_{}x{}_{}.sh".format(percent[j],L,W, iteration), "w")
                    file.write("#!/bin/bash\n")
                    file.write("#$ -N start2\n")
                    file.write("#$ -j y \n")
                    file.write("#$ -l h_rt=12:00:00\n")
                    file.write("#$ -p orucky\n")
                    file.write("#$ -pe mpi_16_tasks_per_node 16\n")
                    file.write("julia -p 16 startenergy{}_{}x{}_{}.jl".format(percent[j], L, W, iteration))
            if L == 10:
                for (i,W) in enumerate([50, 100, 200]):
                    file = open("job_startenergy{}_{}x{}_{}.sh".format(percent[j],L,W, iteration), "w")
                    file.write("#!/bin/bash\n")
                    file.write("#$ -N start2\n")
                    file.write("#$ -j y \n")
                    file.write("#$ -l h_rt=12:00:00\n")
                    file.write("#$ -p orucky\n")
                    file.write("#$ -pe mpi_16_tasks_per_node 16\n")
                    file.write("julia -p 16 startenergy{}_{}x{}_{}.jl".format(percent[j], L, W, iteration))
            if L == 20:
                for (i,W) in enumerate([50, 100, 200]):
                    file = open("job_startenergy{}_{}x{}_{}.sh".format(percent[j],L,W, iteration), "w")
                    file.write("#!/bin/bash\n")
                    file.write("#$ -N start2\n")
                    file.write("#$ -j y \n")
                    file.write("#$ -l h_rt=12:00:00\n")
                    file.write("#$ -p orucky\n")
                    file.write("#$ -pe mpi_16_tasks_per_node 16\n")
                    file.write("julia -p 16 startenergy{}_{}x{}_{}.jl".format(percent[j], L, W, iteration))
                    
