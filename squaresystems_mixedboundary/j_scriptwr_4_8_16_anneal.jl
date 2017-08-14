Concentration = [0, 25, 50, 75, 100]
L_sizes = [4, 8, 16]
Systems_per_size = 32
fixed_percent = 0.6
ps = ["[0.25  , 0.    , 0.25  , 0.25  , 0.25  ]",
      "[0.1875, 0.25  , 0.1875, 0.1875, 0.1875]",
      "[0.125 , 0.5   , 0.125 , 0.125 , 0.125 ]",
      "[0.0625, 0.75  , 0.0625, 0.0625, 0.0625]",
      "[0.    , 1.    , 0.    , 0.    , 0.    ]"]



for (i,perc) in enumerate(Concentration)
  for L in L_sizes
    fixed_left = collect(Int, 1:round(L * fixed_percent, RoundUp))
    for right in 0:L-1
      (L==4  && !(right == 1 || right == 2 )) && continue
      (L==8  && !(right == 4 || right == 5 )) && continue
      (L==16 && !(right == 9 || right == 10)) && continue
      fixed_right = collect(Int, right:right-1+round(L * fixed_percent, RoundUp))
      fixed_right = mod.(fixed_right, L)+1
      for j in 1:(Systems_per_size)
        open("a$(perc)_$(L)x$(L)_right$(right)_60fixed_$(j).jl", "w") do f
          write(f, "push!(LOAD_PATH, pwd())\n")
	  write(f,
"using annealing
using Graphslib
using JLD\n")
	  write(f, "gates = [idid, toff, swsw, swid, idsw]\n")
	  write(f, "file = load(\"$(perc)_$(L)x$(L)_right$(right)_60fixed_$(j).jld\")\n")
          write(f, "ds = file[\"ds\"]\n")
	  write(f, "solution = file[\"solution\"]\n")
	  write(f, "gm = file[\"gm\"]\n")
	  write(f, "nfix = length(file[\"fixed_right\"])\n")
	  write(f, "fixed_right =  hcat(file[\"fixed_right\"], solution[end - ds[1] + file[\"fixed_right\"]])\n")
	  write(f, "fixed_left  =  hcat(file[\"fixed_left\"], solution[file[\"fixed_left\"]])\n")
	  write(f, "g = vertex_lattice(ds, gates, 3, gm=gm)\n")
	  write(f, "init_annealing!(g, ds, fixed_left, fixed_right)\n")
	  write(f, "@time state_counts, energy = anneal!(g,ds,mcsteps=2^8, majority_steps=150, majority_cutoff=0.75, totalsteps =100, calc_energy=true)\n")
	  write(f, "save(\"$(perc)_$(L)x$(L)_right$(right)_60fixed_$(j).jld\", \"ds\", ds, \"solution\", solution, \"fixed_right\", fixed_right, \"fixed_left\", fixed_left, \"state_counts\", state_counts, \"energy\", energy, \"gm\", gm, \"fixed\", g.vs[\"fixed\"])\n")
        end
      end
    end
  end
end
