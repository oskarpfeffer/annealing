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
      open("$(perc)_$(L)x$(L)_right$(right)_60fixed.jl", "w") do f
        write(f, "push!(LOAD_PATH, pwd())\n")
	write(f, "using annealing\n")
	write(f, "using JLD\n")
	write(f, "gates = [idid, toff, swsw, swid, idsw]\n")
        write(f, "ds = [$(L), $(L)]\n")
	write(f, "ps = $(ps[i])\n")
	write(f, "fixed_left = $(fixed_left)\n")
	write(f, "fixed_right = $(fixed_right)\n")
	write(f, "for system in 1:$(Systems_per_size)\n")
	write(f, "  left, right, gm, solution = unique_solution(ds, fixed_left, fixed_right, solution=true, gates=gates, ps=ps, trials=1000)\n")
	write(f, "  save(\"$(perc)_$(L)x$(L)_right$(right)_60fixed_\$(system).jld\", \"fixed_left\", fixed_left, \"fixed_right\", fixed_right, \"solution\", solution, \"gm\", gm, \"ds\", ds)\n")
	write(f, "end")
      end
    end
  end
end

