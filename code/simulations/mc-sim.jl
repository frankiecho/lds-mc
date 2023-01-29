using Revise, DataFrames, Glob, Random, Pipe, ProgressMeter, ThreadsX, Distributed
addprocs(36);
Random.seed!(123456)
@everywhere begin
    using Random
    S = 1:100;
    include("../functions/type-defs.jl")
    include("../functions/optim-functions.jl");
    include("../functions/expected-utility-functions.jl")
    include("../functions/sim-landscape-functions.jl")
    α = 0:0.1:50
    λ = 0:0.02:1
    budget = 100;
end
@everywhere function fcn_mc_sim(i)
    L = fcn_generate_landscape(yy = 10);
    ev_soln = fcn_optim_ev(-L.R; budget = budget);
    ev_rv_soln = fcn_optim_ev(-L.RV; budget = budget);
    ev_ef = EfficiencyFrontier(ev_soln, L.R' * ev_soln, [0], fcn_optim_ev);
    ev_rv_ef = EfficiencyFrontier(ev_rv_soln, L.R' * ev_rv_soln, [0], fcn_optim_ev);
    cvar_ef = fcn_map_ef(L.R, fcn_optim_cvar, budget, λ)
    mstd_ef = fcn_map_ef(L.R, fcn_optim_mstd, budget, λ)
    ef = Result(ev_ef, ev_rv_ef, cvar_ef, mstd_ef)

    uf_type = "CRRA"
    location, scale = fcn_get_location_scale(L.R)
    utility_functions = map(a -> UtilityFunction(uf_type, a, location, scale), α);
    ev_rv_ce = map(u -> fcn_evaluate_ef(ev_rv_ef, u), utility_functions);
    ev_ce = map(u -> fcn_evaluate_ef(ev_ef, u), utility_functions);
    cvar_ce = map(u -> fcn_evaluate_ef(cvar_ef, u), utility_functions)
    mstd_ce = map(u -> fcn_evaluate_ef(mstd_ef, u), utility_functions)
    ce = Result(ev_ce, ev_rv_ce, cvar_ce, mstd_ce);

    summary_func = maximum
    ev_ce_max = map(summary_func, ev_ce);
    ev_rv_ce_max = map(summary_func, ev_rv_ce);
    cvar_ce_max = map(summary_func, cvar_ce);
    mstd_ce_max = map(summary_func, mstd_ce);
    ce_max = Result(ev_ce_max, ev_rv_ce_max, cvar_ce_max, mstd_ce_max)

    return MCResult(ef, ce, ce_max)
end

result = pmap(fcn_mc_sim, S);

ev = @pipe [result[i].ce_max.ev for i=S] |> mapreduce(permutedims, vcat, _);
ev_rv = @pipe [result[i].ce_max.ev_rv for i=S] |> mapreduce(permutedims, vcat, _);
cvar = @pipe [result[i].ce_max.cvar for i=S] |> mapreduce(permutedims, vcat, _);
mstd = @pipe [result[i].ce_max.mstd for i=S] |> mapreduce(permutedims, vcat, _);

result_array = @pipe [ev ./ ev_rv, cvar./ ev_rv, mstd./ ev_rv, cvar ./ mstd] |> map(w -> w.-1,_);
result_mean = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 50),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
result_lb = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 5),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
result_ub = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 95),x,dims=1),_) |> mapreduce(permutedims, hcat, _);

result_df = DataFrame(vcat(result_mean, result_lb, result_ub), :auto);
rename!(result_df, [:ev,:cvar, :mstd, :cvar_mstd])

result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)));
result_df.alpha = repeat(α, 3)
CSV.write("output/ce_df.csv", result_df)