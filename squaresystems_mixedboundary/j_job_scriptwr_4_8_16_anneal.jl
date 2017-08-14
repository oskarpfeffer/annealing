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
        open("a$(perc)_$(L)x$(L)_right$(right)_60fixed_$(j).sh", "w") do f
          write(f,
"#!/bin/bash
#\$ -N j$(perc)_$(L)_$(j)
#\$ -j y 
#\$ -l h_rt=12:00:00
#\$ -p orucky
julia a$(perc)_$(L)x$(L)_right$(right)_60fixed_$(j).jl")
        end
      end
    end
  end
end
