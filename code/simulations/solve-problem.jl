includet("../functions/optim-functions.jl");
includet("generate-landscape.jl");

using StatsPlots;
plotlyjs()

budget = 200;
cvar_soln = fcn_optim_cvar(-Wps; budget = budget, λ=1);
mv_soln = fcn_optim_mv(-Wps; budget = budget, λ=1);
ev_soln = fcn_optim_ev(-Wps; budget = budget);
soln = hcat(ev_soln, cvar_soln, mv_soln);
returns = Wps' * soln;
fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
plot_vecs = hcat(mean(Wps, dims=2), mean(Wps, dims = 2), pw, sum(SS, dims=2), ev_soln, cvar_soln);
plots = [fcn_plot_heatmap(plot_vecs[:,i]) for i in axes(plot_vecs,2)];

hist_plot = histogram(returns, layout = (3,1), xlims = (minimum(returns), maximum(returns)), nbins = 30);
heatmap_plot = plot(plots...; aspect_ratio=:equal);

hist_plot