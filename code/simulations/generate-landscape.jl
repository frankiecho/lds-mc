using Revise

includet("../functions/sim-landscape-functions.jl")

Random.seed!(1234);
dims = (50,50);
S = 100;
W, P = fcn_spatial_weights(dims; distance = true, bandwidth = 5);
X=1 .+rand(size(W,1),1);
Wp = fcn_spatial_AR(W, S; X = X, ρ=0.9, σ=5);
Wμ = mean(Wp,dims=2)
# x_hms = heatmap(reshape(X, dims));hms = [heatmap(reshape(Wp[:,s], dims)) for s=1:4];
# hms_mean = heatmap(reshape(mean(Wp,dims=2), dims));
# plot(hms...;aspect_ratio = :equal, grid=false);
# plot(hms_mean)

n_shocks = Poisson(1);
shock_size = Poisson(20);
pw = (Wμ.-minimum(Wμ) + 2*rand(size(Wμ,1)))  |> vec; # Probability of shock
SS = fcn_spatial_shock(W, S, n_shocks, shock_size, p=pw);
Wps = Wp;
Wps[SS.==1] .= 0;

fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
fcn_plot_heatmap(mean(Wps, dims=2))