include("../../code/functions/type-defs.jl")
include("../../code/functions/mc-sim-functions.jl")
α = 0:0.1:100
test_param = LandscapeParameters((10,10), 0, 0.5, 0.5, 0.0001, 0.9)
result = map((ii) -> fcn_mc_sim(ii, test_param; λ = 0:0.05:1, α=α), 1:1);
fcn_write_contiguity(result, "param_search", true, 0:0.1:100)