using Distributed, JLD2

addprocs(35);
@everywhere using Gurobi
using Gurobi

#GRB_ENV = Gurobi.Env(); sendto(workers(),GRB_ENV=GRB_ENV);
@everywhere begin
    using Revise, DataFrames, Glob, Random, Pipe, ProgressMeter, CSV
    Random.seed!(123456)
    using Random
    Q = 1:50;
    #cd("/Users/frankiecho/Documents/Github/lds-mc-julia/code/simulations/")
    include("../../code/functions/mc-sim-functions.jl")
end

result = pmap(i -> fcn_mc_sim(i), Q);
fcn_write_result(result, "baseline");
fcn_write_shock_exposure(result, "baseline");
fcn_write_contiguity(result, "baseline");

#result_no_risk = pmap(i -> fcn_mc_sim(i, 0), Q);
#fcn_write_result(result_no_risk, "norisk"); 

#fcn_save_results([result]);