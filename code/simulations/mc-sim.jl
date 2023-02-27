using Distributed, JLD2

addprocs(5);
@everywhere using Gurobi
using Gurobi

#GRB_ENV = Gurobi.Env(); sendto(workers(),GRB_ENV=GRB_ENV);
@everywhere begin
using Revise, DataFrames, Glob, Random, Pipe, ProgressMeter, CSV
Random.seed!(123456)
using Random
Q = 1:5;
#cd("/Users/frankiecho/Documents/Github/lds-mc-julia/code/simulations/")
include("../../code/functions/mc-sim-functions.jl")
end

result = pmap(i -> fcn_mc_sim(i), Q);
fcn_write_result(result, "baseline");

result_no_risk = pmap(i -> fcn_mc_sim(i, 0), Q);
fcn_write_result(result_no_risk, "norisk"); 

#jldsave("../../output/result.jld2"; result, result_no_risk)