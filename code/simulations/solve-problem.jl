using Revise

includet("../functions/optim-functions.jl");
includet("../functions/expected-utility-functions.jl");
includet("generate-landscape.jl");

using StatsPlots;
plotlyjs()
theme(:juno)

budget = 100;
map_ef = (f,lambda_vec,R) -> pmap(l->f(-R; budget = budget, λ=l),lambda_vec) |> e->mapreduce(permutedims, vcat, e) |> transpose;
plot_ef = (returns, risk_func) -> hcat(mean(returns, dims = 1)', risk_func(returns)) |> m->plot!(m[:,2],m[:,1], seriestype=:scatter);
plot_dens = (returns) -> density(returns, labels = string.(lambda_vec'));

lambda_vec = 0:0.1:1;
mstd_soln = map_ef(fcn_optim_mstd, lambda_vec, Wps);
#mv_soln = map_ef(fcn_optim_mv, lambda_vec, Wps);
cvar_soln = map_ef(fcn_optim_cvar, lambda_vec, Wps);
ev_soln = fcn_optim_ev(-Wps; budget = budget);
risk_threshold = percentile(Wps' * ev_soln, 10);
scaling_factor = percentile(Wps' * ev_soln, 90)-percentile(Wps' * ev_soln, 10);
EU_scale = (w)->(w.-risk_threshold)./scaling_factor;
EU_scale_inv = (w) -> (w .* scaling_factor) .+ risk_threshold;

mstd_returns = Wps' * mstd_soln;
cvar_returns = Wps' * cvar_soln;
plot()
plot_dens(mstd_returns)
plot_dens(cvar_returns)

find_max_α = (α, w) -> maximum(EU(EU_scale(w),v->CARA(v,α))) |> w->CARA_inv(w,α) |> EU_scale_inv
find_max_γ = (γ, w) -> maximum(EU(EU_scale(w).+1,v->CRRA(v,γ))) |> w->CRRA_inv(w,γ) |> EU_scale_inv
find_max_u_pl = (r, w) -> maximum(mean(piecewise_linear(w, risk_threshold, r, 1), dims = 1))

α_vec = 1:0.1:30
mstd_ce_max = map(α -> find_max_α(α, mstd_returns), α_vec);
cvar_ce_max = map(α -> find_max_α(α, cvar_returns), α_vec);

mstd_ce_max_crra = map(γ -> find_max_γ(γ, mstd_returns), α_vec);
cvar_ce_max_crra = map(γ -> find_max_γ(γ, cvar_returns), α_vec);

mstd_pl_max = map(α -> find_max_u_pl(α, mstd_returns), α_vec);
cvar_pl_max = map(α -> find_max_u_pl(α, cvar_returns), α_vec);

plot(α_vec, [mstd_pl_max, cvar_pl_max], labels = ["M-SD" "M-CVaR"])
plot(α_vec, [mstd_ce_max, cvar_ce_max], labels = ["M-SD" "M-CVaR"])
plot(α_vec, [mstd_ce_max_crra, cvar_ce_max_crra], labels = ["M-SD" "M-CVaR"])


α = 4;
fcn_ce = w->EU(w|>EU_scale,v->CARA(v,α));
fcn_ce_pl = w->EU(w, v->piecewise_linear(v,risk_threshold,10,0))
plot()
plot_ef(Wps' * mstd_soln, fcn_ce)
plot_ef(Wps' * cvar_soln, fcn_ce)

fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
plot_vecs = hcat(mean(Wp, dims=2), mean(Wps, dims = 2), pw, sum(SS, dims=2), ev_soln, cvar_soln[:,11]);
plots = [fcn_plot_heatmap(plot_vecs[:,i]) for i in axes(plot_vecs,2)];
heatmap_plot = plot(plots...; aspect_ratio=:equal);
heatmap_plot

returns = Wps' * mstd_soln;
hist_plot = histogram(returns, layout = (3,1), xlims = (minimum(returns), maximum(returns)), nbins = 30);


hist_plot

