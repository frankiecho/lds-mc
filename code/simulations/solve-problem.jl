using Revise

includet("../functions/optim-functions.jl");
includet("generate-landscape.jl");

using StatsPlots;
plotlyjs()

budget = 200;
map_ef = (f,lambda_vec,R) -> pmap(l->f(-R; budget = budget, λ=l),lambda_vec) |> e->mapreduce(permutedims, vcat, e) |> transpose;
plot_ef = (returns) -> [mean(returns, dims = 1); std(returns, dims = 1)]' |> m->plot(m[:,2],m[:,1], seriestype=:scatter);

lambda_vec = 0:0.1:1;
mstd_soln = map_ef(fcn_optim_mstd, lambda_vec, Wps);
mv_soln = map_ef(fcn_optim_mv, lambda_vec, Wps);
cvar_soln = map_ef(fcn_optim_cvar, lambda_vec, Wps);
plot_ef(Wps' * mstd_soln)

plot(x->CARA(x,1), -3, 5)

cvar_soln = pmap(l->fcn_optim_cvar(-Wps; budget = budget, λ=l),1) |> e->mapreduce(permutedims, vcat, e);
mv_soln = fcn_optim_mv(-Wps; budget = budget, λ=0);
ev_soln = fcn_optim_ev(-Wps; budget = budget);
soln = hcat(ev_soln, cvar_soln', mv_soln, mstd_soln);
returns = Wps' * soln;
fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
plot_vecs = hcat(mean(Wps, dims=2), mean(Wps, dims = 2), pw, sum(SS, dims=2), ev_soln, cvar_soln);
plots = [fcn_plot_heatmap(plot_vecs[:,i]) for i in axes(plot_vecs,2)];

hist_plot = histogram(returns, layout = (3,1), xlims = (minimum(returns), maximum(returns)), nbins = 30);
heatmap_plot = plot(plots...; aspect_ratio=:equal);

hist_plot

