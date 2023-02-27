home_dir = "/users/frankiecho/Documents/Github/lds-mc-julia"
include("$(home_dir)/code/functions/type-defs.jl")
include("$(home_dir)/code/functions/optim-functions.jl");
include("$(home_dir)/code/functions/expected-utility-functions.jl")
include("$(home_dir)/code/functions/sim-landscape-functions.jl")

α = 0:0.1:50
λ = 0:0.1:1
budget = 100;


using Pipe

function fcn_mc_sim(i, n_shocks=1000)
    L = fcn_generate_landscape((40,40), n_shocks=n_shocks)
    #L = fcn_generate_landscape((10,10), yy = 10, n_shocks=n_shocks, shock_size=Poisson(1));
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

    println("Completed run $(i)")
    return MCResult(ef, ce, ce_max, L)
end

function fcn_write_result(result::AbstractArray, suffix = "")
    ev = @pipe [result[q].ce_max.ev for q=Q] |> mapreduce(permutedims, vcat, _);
    ev_rv = @pipe [result[q].ce_max.ev_rv for q=Q] |> mapreduce(permutedims, vcat, _);
    cvar = @pipe [result[q].ce_max.cvar for q=Q] |> mapreduce(permutedims, vcat, _);
    mstd = @pipe [result[q].ce_max.mstd for q=Q] |> mapreduce(permutedims, vcat, _);

    result_array = @pipe [ev ./ ev_rv, cvar./ ev_rv, mstd./ ev_rv, cvar ./ mstd] |> map(w -> w.-1,_);
    result_mean = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 50),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
    result_lb = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 5),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
    result_ub = @pipe result_array |> map(x->mapslices(xx->percentile(xx, 95),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
    result_min = @pipe result_array |> map(x->mapslices(xx->minimum(xx),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
    result_max = @pipe result_array |> map(x->mapslices(xx->maximum(xx),x,dims=1),_) |> mapreduce(permutedims, hcat, _);
    result_df = DataFrame(vcat(result_mean, result_lb, result_ub, result_min, result_max), :auto);
    rename!(result_df, [:ev,:cvar, :mstd, :cvar_mstd])

    result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)), repeat(["min"], length(α)), repeat(["max"], length(α)));
    result_df.alpha = repeat(α, 5)
    CSV.write("$(home_dir)/output/ce_df_$(suffix)_$(Q[end]).csv", result_df)
end
