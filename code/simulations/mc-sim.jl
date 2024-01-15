using Distributed, JLD2

addprocs(22);
@everywhere using Gurobi
using Gurobi

#GRB_ENV = Gurobi.Env(); sendto(workers(),GRB_ENV=GRB_ENV);
@everywhere begin
    using Revise, DataFrames, Glob, Random, Pipe, ProgressMeter, CSV
    Random.seed!(123456)
    using Random
    
    #cd("/Users/frankiecho/Documents/Github/lds-mc-julia/code/simulations/")
    include("../../code/simulations/specify-params.jl")
    include("../../code/functions/mc-sim-functions.jl")
    Q = 1:length(param_vec);
    
    function fcn_mc_sim_i(i)
        α = 0:0.1:100
        result = map((ii) -> fcn_mc_sim(ii, param_vec[i]; λ = 0:0.05:1, α=α), 1:nsims);
        fcn_write_result(result, "param_search_$(i)"; α=α);
        fcn_write_downside(result, "param_search_$(i)", 0:0.1:100);
        fcn_write_contiguity(result, "param_search_$(i)",  true, 0:0.1:100);
        println("Completed run $(i)")
        return
    end
end


result = pmap(i -> fcn_mc_sim_i(i), Q);

#result_no_risk = pmap(i -> fcn_mc_sim(i, 0), Q);
#fcn_write_result(result_no_risk, "norisk"); 

#fcn_save_results([result]);