using Revise

includet("../functions/sim-landscape-functions.jl")

Random.seed!(1234);
dims = (50,50);
S = 2501;
n_shocks = Poisson(1);
shock_size = Poisson(20);

W, P = fcn_spatial_weights(dims; distance = true, bandwidth = 5);
X    = 20 .+rand(size(W,1),1);
Wp   = fcn_spatial_AR(W, S; X = X, ρ=0.9, σ=5);
Wp[Wp .< 0] .= 0;
Wμ   = mean(Wp,dims=2)

pw   = (Wμ.-minimum(Wμ) + 2*rand(size(Wμ,1))) |> vec; # Probability of shock
SS   = fcn_spatial_shock(W, S, n_shocks, shock_size, p=pw);
Wps  = Wp;
Wps[SS.==1] .= 0;

#fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
#fcn_plot_heatmap(mean(Wps, dims=2))