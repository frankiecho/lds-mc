using JLD2

home_dir = "d:/Github/lds-mc"
include("$(home_dir)/code/functions/mc-sim-functions.jl")
result = load("$(home_dir)/output/result.jld2");

result = result["result_array"][1];
fcn_write_shock_exposure(result, "baseline", 30:2:80, 1:100);
fcn_write_contiguity(result, "baseline", 1:100);