push!(LOAD_PATH, pwd())
using annealing 
using Graphslib
using JLD
gates = [idid, toff, swsw, swid, idsw]
file = load("75_16x16_right10_60fixed_32.jld")
ds = file["ds"]
solution = file["solution"]
gm = file["gm"]
nfix = length(file["fixed_right"])
fixed_right =  hcat(file["fixed_right"], solution[end - ds[1] + file["fixed_right"]])
fixed_left  =  hcat(file["fixed_left"], solution[file["fixed_left"]])
g = vertex_lattice(ds, gates, 3, gm=gm)
init_annealing!(g, ds, fixed_left, fixed_right)
@time state_counts, energy = anneal!(g,ds,mcsteps=2^8, majority_steps=150, majority_cutoff=0.75, totalsteps =100, calc_energy=true)

save("a75_16x16_right10_60fixed_32.jld", "ds", ds, "solution", solution, "fixed_right", fixed_right, "fixed_left", fixed_left, "state_counts", state_counts, "energy", energy, "gm", gm, "fixed", g.vs["fixed"])
