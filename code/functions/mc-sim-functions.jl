#home_dir = "/users/frankiecho/Documents/Github/lds-mc"
home_dir = "d:/Github/lds-mc"
include("$(home_dir)/code/functions/type-defs.jl")
include("$(home_dir)/code/functions/optim-functions.jl");
include("$(home_dir)/code/functions/expected-utility-functions.jl")
include("$(home_dir)/code/functions/sim-landscape-functions.jl")

using CSV
using Pipe

function fcn_mc_sim(i, params::LandscapeParameters=LandscapeParameters((40, 40), 10, 0.7, 5, 0.0001, 100, 0.99, (false, true), 0:0.1:1); α=0:0.1:50)
    w0 = params.budget * params.yy + 1
    budget = params.budget
    λ = params.λ
    L = fcn_generate_landscape(params.dims; yy=params.yy, ρ=params.ρ, σ=params.σ, η=params.η, shocks=params.shocks)
    ev_soln = fcn_optim_ev(-L.R; budget=budget)
    ev_rv_soln = fcn_optim_ev(-L.RV; budget=budget)
    ev_ef = EfficiencyFrontier(ev_soln, w0 .+ L.R' * ev_soln, [0], fcn_optim_ev)
    ev_rv_ef = EfficiencyFrontier(ev_rv_soln, w0 .+ L.R' * ev_rv_soln, [0], fcn_optim_ev)

    # Find optimal lambda for each problem
    #λ_cvar = map(u -> fcn_max_lambda(u, L.R, fcn_optim_cvar, params.budget, params.β, w0), utility_functions)
    #λ_mstd = map(u -> fcn_max_lambda(u, L.R, fcn_optim_mstd, params.budget, 0, w0), utility_functions)

    cvar_ef = fcn_map_ef(L.R, fcn_optim_cvar, budget, λ, params.β, w0)
    mstd_ef = fcn_map_ef(L.R, fcn_optim_mstd, budget, λ, 0, w0)
    ef = Result(ev_ef, ev_rv_ef, cvar_ef, mstd_ef)

    uf_type = "CRRA"
    #location, scale = fcn_get_location_scale(L.R)
    w = mean(L.R' * ev_soln)
    utility_functions = map(a -> UtilityFunction(uf_type, a, 0, w), α)
    ev_rv_ce = map(u -> fcn_evaluate_ef(ev_rv_ef, u), utility_functions)
    ev_ce = map(u -> fcn_evaluate_ef(ev_ef, u), utility_functions)
    cvar_ce = map(u -> fcn_evaluate_ef(cvar_ef, u), utility_functions)
    mstd_ce = map(u -> fcn_evaluate_ef(mstd_ef, u), utility_functions)
    ce = Result(ev_ce, ev_rv_ce, cvar_ce, mstd_ce)

    summary_func = maximum
    ev_ce_max = map(summary_func, ev_ce)
    ev_rv_ce_max = map(summary_func, ev_rv_ce)
    cvar_ce_max = map(summary_func, cvar_ce)
    mstd_ce_max = map(summary_func, mstd_ce)
    ce_max = Result(ev_ce_max, ev_rv_ce_max, cvar_ce_max, mstd_ce_max)


    return MCResult(ef, ce, ce_max, L)
end

function fcn_write_result(result::AbstractArray, suffix=""; α=0:0.1:50)
    Q = 1:length(result)
    ev = @pipe [result[q].ce_max.ev for q = Q] |> mapreduce(permutedims, vcat, _)
    ev_rv = @pipe [result[q].ce_max.ev_rv for q = Q] |> mapreduce(permutedims, vcat, _)
    cvar = @pipe [result[q].ce_max.cvar for q = Q] |> mapreduce(permutedims, vcat, _)
    mstd = @pipe [result[q].ce_max.mstd for q = Q] |> mapreduce(permutedims, vcat, _)

    cols = [:ev, :cvar, :mstd, :cvar_ev, :mstd_ev, :cvar_mstd]
    result_array = @pipe [ev ./ ev_rv, cvar ./ ev_rv, mstd ./ ev_rv, cvar ./ ev, mstd ./ ev, cvar ./ mstd] |> map(w -> w .- 1, _)
    result_mean = @pipe result_array |> map(x -> mapslices(xx -> percentile(xx, 50), x, dims=1), _) |> mapreduce(permutedims, hcat, _)
    result_lb = @pipe result_array |> map(x -> mapslices(xx -> percentile(xx, 5), x, dims=1), _) |> mapreduce(permutedims, hcat, _)
    result_ub = @pipe result_array |> map(x -> mapslices(xx -> percentile(xx, 95), x, dims=1), _) |> mapreduce(permutedims, hcat, _)
    result_min = @pipe result_array |> map(x -> mapslices(xx -> minimum(xx), x, dims=1), _) |> mapreduce(permutedims, hcat, _)
    result_max = @pipe result_array |> map(x -> mapslices(xx -> maximum(xx), x, dims=1), _) |> mapreduce(permutedims, hcat, _)
    result_df = DataFrame(vcat(result_mean, result_lb, result_ub, result_min, result_max), :auto)
    rename!(result_df, cols)
    result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)), repeat(["min"], length(α)), repeat(["max"], length(α)))
    result_df.alpha = repeat(α, 5)
    CSV.write("$(home_dir)/output/ce_df_$(suffix)_$(Q[end]).csv", result_df)
end

function fcn_write_lambda(result::AbstractArray, suffix=""; α=0:0.1:50, λ=0:0.1:1)
    Q = 1:length(result)
    cvar_lambda = hcat([map(e -> λ[argmax(e)[2]], res.ce.cvar) for res in result]...)
    mstd_lambda = hcat([map(e -> λ[argmax(e)[2]], res.ce.mstd) for res in result]...)
    cvar_lambda_mean = mapslices(xx -> percentile(xx, 50), cvar_lambda, dims=2)
    cvar_lambda_lb = mapslices(xx -> percentile(xx, 5), cvar_lambda, dims=2)
    cvar_lambda_ub = mapslices(xx -> percentile(xx, 95), cvar_lambda, dims=2)

    mstd_lambda_mean = mapslices(xx -> percentile(xx, 50), mstd_lambda, dims=2)
    mstd_lambda_lb = mapslices(xx -> percentile(xx, 5), mstd_lambda, dims=2)
    mstd_lambda_ub = mapslices(xx -> percentile(xx, 95), mstd_lambda, dims=2)

    cols = [:cvar_mean, :cvar_lb, :cvar_ub, :mstd_mean, :mstd_lb, :mstd_ub]
    result_df = DataFrame(hcat(cvar_lambda_mean, cvar_lambda_lb, cvar_lambda_ub, mstd_lambda_mean, mstd_lambda_lb, mstd_lambda_ub), :auto)
    rename!(result_df, cols)
    result_df.alpha = α
    CSV.write("$(home_dir)/output/lambda_df_$(suffix)_$(Q[end]).csv", result_df)
end

function fcn_write_shock_exposure(result::AbstractArray, suffix="", threshold=50:2:100, Q=1:5)
    for t = threshold
        likelihood_shock_mstd = zeros(length(α), length(Q))
        likelihood_shock_cvar = zeros(length(α), length(Q))
        for q = Q
            max_ce_id_mstd = map(i -> findmax(i)[2][2], result[q].ce.mstd)
            max_ce_id_cvar = map(i -> findmax(i)[2][2], result[q].ce.cvar)

            ev_rv_likelihood = mean(((result[q].L.R .< 1e-5)' * result[q].ef.ev_rv.solutions) .> t)
            mstd_likelihood = (mean(((result[q].L.R .< 1e-5)' * result[q].ef.mstd.solutions) .> t, dims=2))
            cvar_likelihood = (mean(((result[q].L.R .< 1e-5)' * result[q].ef.mstd.solutions) .> t, dims=2))

            for (i, x) = enumerate(max_ce_id_mstd)
                likelihood_shock_mstd[i, q] = mstd_likelihood[x]
            end

            for (i, x) = enumerate(max_ce_id_cvar)
                likelihood_shock_cvar[i, q] = cvar_likelihood[x]
            end
        end

        run_median_mstd = percentile.(eachrow(likelihood_shock_mstd), 50)
        run_ub_mstd = percentile.(eachrow(likelihood_shock_mstd), 95)
        run_lb_mstd = percentile.(eachrow(likelihood_shock_mstd), 5)

        run_median_cvar = percentile.(eachrow(likelihood_shock_cvar), 50)
        run_ub_cvar = percentile.(eachrow(likelihood_shock_cvar), 95)
        run_lb_cvar = percentile.(eachrow(likelihood_shock_cvar), 5)

        result_df = DataFrame(hcat(vcat(run_median_mstd, run_lb_mstd, run_ub_mstd), vcat(run_median_cvar, run_lb_cvar, run_ub_cvar)), :auto)
        rename!(result_df, [:mstd, :cvar])
        result_df.alpha = repeat(α, 3)
        result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)))
        CSV.write("$(home_dir)/output/likelihood_$(suffix)_$(Q[end])_$(t).csv", result_df)
    end
end

function fcn_write_contiguity(result::AbstractArray, suffix="", distance=false, α=0:0.1:50)
    Q = 1:length(result)
    contiguity_mstd = zeros(length(α), length(Q))
    contiguity_cvar = zeros(length(α), length(Q))
    q = 1
    dims = result[q].L.dims
    W, P, b, D = fcn_spatial_weights(dims; boundary=1, distance=true, bandwidth=sqrt(sum(dims .^ 2)), α=2)
    D = D[b.==0, b.==0]

    for q = Q
        max_ce_id_mstd = map(i -> findmax(i)[2][2], result[q].ce.mstd)
        max_ce_id_cvar = map(i -> findmax(i)[2][2], result[q].ce.cvar)

        if (distance)
            contiguity_ev_lambda = result[q].ef.ev.solutions' * D * result[q].ef.ev.solutions
            contiguity_mstd_lambda = (map((i) -> result[q].ef.mstd.solutions[:, i]' * D * result[q].ef.mstd.solutions[:, i], 1:size(result[q].ef.mstd.solutions, 2)) .- contiguity_ev_lambda) ./ contiguity_ev_lambda
            contiguity_cvar_lambda = (map((i) -> result[q].ef.cvar.solutions[:, i]' * D * result[q].ef.cvar.solutions[:, i], 1:size(result[q].ef.cvar.solutions, 2)) .- contiguity_ev_lambda) ./ contiguity_ev_lambda
        else
            contiguity_ev_lambda = result[q].ef.ev.solutions' * result[q].L.W * result[q].ef.ev.solutions
            contiguity_mstd_lambda = (map((i) -> result[q].ef.mstd.solutions[:, i]' * result[q].L.W * result[q].ef.mstd.solutions[:, i], 1:size(result[q].ef.mstd.solutions, 2)) .- contiguity_ev_lambda) ./ contiguity_ev_lambda
            contiguity_cvar_lambda = (map((i) -> result[q].ef.cvar.solutions[:, i]' * result[q].L.W * result[q].ef.cvar.solutions[:, i], 1:size(result[q].ef.cvar.solutions, 2)) .- contiguity_ev_lambda) ./ contiguity_ev_lambda
        end
        for (i, x) = enumerate(max_ce_id_mstd)
            contiguity_mstd[i, q] = contiguity_mstd_lambda[x]
        end

        for (i, x) = enumerate(max_ce_id_cvar)
            contiguity_cvar[i, q] = contiguity_cvar_lambda[x]
        end
    end
    run_median_mstd = percentile.(eachrow(contiguity_mstd), 50)
    run_ub_mstd = percentile.(eachrow(contiguity_mstd), 95)
    run_lb_mstd = percentile.(eachrow(contiguity_mstd), 5)

    run_median_cvar = percentile.(eachrow(contiguity_cvar), 50)
    run_ub_cvar = percentile.(eachrow(contiguity_cvar), 95)
    run_lb_cvar = percentile.(eachrow(contiguity_cvar), 5)

    result_df = DataFrame(hcat(vcat(run_median_mstd, run_lb_mstd, run_ub_mstd), vcat(run_median_cvar, run_lb_cvar, run_ub_cvar)), :auto)
    rename!(result_df, [:mstd, :cvar])
    result_df.alpha = repeat(α, 3)
    result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)))
    if (distance)
        CSV.write("$(home_dir)/output/distance_$(suffix)_$(Q[end]).csv", result_df)
    else
        CSV.write("$(home_dir)/output/contiguity_$(suffix)_$(Q[end]).csv", result_df)
    end
end


function fcn_write_downside(result::AbstractArray, suffix="", α=0:0.1:50, pct=(5, 95))
    Q = 1:length(result)
    mstd_downside_crra = zeros(length(α), length(Q))
    cvar_downside_crra = zeros(length(α), length(Q))
    mstd_upside_crra = zeros(length(α), length(Q))
    cvar_upside_crra = zeros(length(α), length(Q))
    for q = Q
        max_ce_id_mstd = map(i -> findmax(i)[2][2], result[q].ce.mstd)
        max_ce_id_cvar = map(i -> findmax(i)[2][2], result[q].ce.cvar)
        lb = percentile(result[q].ef.ev.returns, pct[1])
        ub = percentile(result[q].ef.ev.returns, pct[2])

        cvar_downside = mean(result[q].ef.cvar.returns .<= lb, dims=1)'
        mstd_downside = mean(result[q].ef.mstd.returns .<= lb, dims=1)'
        cvar_upside = mean(result[q].ef.cvar.returns .>= ub, dims=1)'
        mstd_upside = mean(result[q].ef.mstd.returns .>= ub, dims=1)'

        for (i, x) = enumerate(max_ce_id_mstd)
            mstd_downside_crra[i, q] = mstd_downside[x]
            mstd_upside_crra[i, q] = mstd_upside[x]
        end

        for (i, x) = enumerate(max_ce_id_cvar)
            cvar_downside_crra[i, q] = cvar_downside[x]
            cvar_upside_crra[i, q] = cvar_upside[x]
        end
    end

    run_median_downside_mstd = percentile.(eachrow(mstd_downside_crra), 50)
    run_ub_downside_mstd = percentile.(eachrow(mstd_downside_crra), 95)
    run_lb_downside_mstd = percentile.(eachrow(mstd_downside_crra), 5)

    run_median_upside_mstd = percentile.(eachrow(mstd_upside_crra), 50)
    run_ub_upside_mstd = percentile.(eachrow(mstd_upside_crra), 95)
    run_lb_upside_mstd = percentile.(eachrow(mstd_upside_crra), 5)

    run_median_downside_cvar = percentile.(eachrow(cvar_downside_crra), 50)
    run_ub_downside_cvar = percentile.(eachrow(cvar_downside_crra), 95)
    run_lb_downside_cvar = percentile.(eachrow(cvar_downside_crra), 5)

    run_median_upside_cvar = percentile.(eachrow(cvar_upside_crra), 50)
    run_ub_upside_cvar = percentile.(eachrow(cvar_upside_crra), 95)
    run_lb_upside_cvar = percentile.(eachrow(cvar_upside_crra), 5)

    result_df = DataFrame(hcat(vcat(run_median_downside_mstd, run_lb_downside_mstd, run_ub_downside_mstd), vcat(run_median_downside_cvar, run_lb_downside_cvar, run_ub_downside_cvar)), :auto)
    rename!(result_df, [:mstd, :cvar])
    result_df.alpha = repeat(α, 3)
    result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)))
    CSV.write("$(home_dir)/output/downside_$(suffix)_$(Q[end]).csv", result_df)

    result_df = DataFrame(hcat(vcat(run_median_upside_mstd, run_lb_upside_mstd, run_ub_upside_mstd), vcat(run_median_upside_cvar, run_lb_upside_cvar, run_ub_upside_cvar)), :auto)
    rename!(result_df, [:mstd, :cvar])
    result_df.alpha = repeat(α, 3)
    result_df.var = vcat(repeat(["median"], length(α)), repeat(["lb"], length(α)), repeat(["ub"], length(α)))
    CSV.write("$(home_dir)/output/upside_$(suffix)_$(Q[end]).csv", result_df)
end

function fcn_save_results(result_array)
    jldsave("$(home_dir)/output/result.jld2"; result_array)
end