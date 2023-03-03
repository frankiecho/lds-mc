using Revise
cd("/Users/frankiecho/Documents/Github/lds-mc-julia/code/simulations/")
include("../functions/type-defs.jl")
include("../functions/optim-functions.jl");
include("../functions/expected-utility-functions.jl")
include("../functions/sim-landscape-functions.jl")

Random.seed!(1234);
dims = (40,40);
S = prod(dims) + 1;
n_shocks = 50;
shock_size = Poisson(20);
yy = 200;

W, P, b, D = fcn_spatial_weights(dims; distance = true, bandwidth = 20, boundary = 1, α=2);
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
fcn_plot_heatmap = D -> D |> (x->fcn_reshape_lds(x, dims)) |> heatmap;
s = 10;
fcn_plot_heatmap(Wp[:,s])
fcn_plot_heatmap(SS[:,s])
fcn_plot_heatmap(Wps[:,s])
#fcn_plot_heatmap(Vector(W[:,500]);


using Plots, Distributions
L = fcn_generate_landscape(yy = 10)
histogram(L.nss)
suffstats(L.nss)

SS, NSS, sl = fcn_spatial_shock_gev(W, 500, p = 0.0005, scale = 5);
fcn_plot_heatmap(SS[b .== 0,30])


W, P, b, D = fcn_spatial_weights((40,40); boundary = 1, distance = true, bandwidth = 20,  α=2);
X    = rand(Normal(), size(W,1));
Wp   = fcn_spatial_AR(W, S; X = X, γ=0.9, θ=0.2, σ=1);
p1 = fcn_plot_heatmap(Wp[b .== 0,20])
p2 = fcn_plot_heatmap(Wp[b .== 0,300])
p3 = fcn_plot_heatmap(mean(Wp[b .== 0,:], dims =2))
plot(p1,p2, p3,size=(700,500), layout = (2,2))



result.