push!(LOAD_PATH, pwd())
using annealing
using JLD
gates = [idid, toff, swsw, swid, idsw]
ds = [30,5]
ps=[2.,2.,2.,2.,2.]
decades = 12
totalsteps = 2 
majoritysteps = 256 
fixed_left = [1,2,3]
fixed_right = [1,2,3]
statecounts("compile", 1, ds, gates, ps, 1, 2, fixed_left , fixed_right)
@time statecounts("dat/fixedboundary/30x5_3left3right_12dec_2steps_20toff", decades
, ds, gates, ps, totalsteps, majoritysteps, fixed_left , fixed_right)
