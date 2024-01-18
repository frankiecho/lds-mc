using Revise, DataFrames, Glob, Random, Pipe, ProgressMeter, CSV, Plots
#cd("/Users/frankiecho/Documents/Github/lds-mc")
include("../../code/functions/type-defs.jl")
include("../../code/functions/optim-functions.jl");
include("../../code/functions/expected-utility-functions.jl")
include("../../code/functions/sim-landscape-functions.jl")
include("../../code/functions/mc-sim-functions.jl")
include("../../code/simulations/specify-params.jl")
p = param_vec[1]
p.yy = 10
p.η = 0.01
p.β = 0.99
p.dims = (20, 20)
res = fcn_mc_sim(1, p)