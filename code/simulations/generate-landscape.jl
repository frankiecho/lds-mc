using Revise
includet("../functions/sim-landscape-functions.jl")

Random.seed!(1234);
dims = (40,40);
S = prod(dims) + 1;
n_shocks = Poisson(1);
shock_size = Poisson(20);
yy = 200;

W, P, b = fcn_spatial_weights(dims; distance = true, bandwidth = 20, boundary = 1, α=2);
X    = rand(Normal(), size(W,1));
Wp   = fcn_spatial_AR(W, S; X = X, ρ=0.5, σ=5) .+ yy;
Wp = Wp[b .== 0, :];
Wp[Wp .< 0] .= 0;
Wμ   = mean(Wp,dims=2);

pw   = fcn_spatial_AR(W, 1; X = rand(size(W,1),1), ρ=0.9, σ=5) |> vec; # Probability of shock
pw   = pw[b .== 0,:];
pw[pw .< 0] .= 0;
pw   = pw ./ sum(pw);
SS, nss, sl = fcn_spatial_shock(W[b .== 0, b .== 0], S, n_shocks, shock_size, p=vec(pw));
Wps  = copy(Wp);
Wps[SS .> 1e-5] .= 0;

fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> d->heatmap(d; yflip=true, dpi =600, size=(750,750));
#s = 3;
#fcn_plot_heatmap(Wp[:,s])
#fcn_plot_heatmap(SS[:,s])
#fcn_plot_heatmap(Wps[:,s])

#fcn_plot_heatmap(Vector(W[:,500]))

q, s = findmax(nss);
fcn_plot_heatmap(Wp[:,s]) |> f->savefig(f,"plots/Wp.png");
fcn_plot_heatmap(SS[:,s]) |> f->savefig(f,"plots/SS.png");
fcn_plot_heatmap(Wps[:,s]) |> f->savefig(f,"plots/Wps.png");
=======
fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
s = 10;
fcn_plot_heatmap(Wp[:,s])
fcn_plot_heatmap(SS[:,s])
fcn_plot_heatmap(Wps[:,s])
#fcn_plot_heatmap(Vector(W[:,500]))
>>>>>>> Stashed changes
