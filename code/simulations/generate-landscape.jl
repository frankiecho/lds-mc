using Revise

includet("../functions/sim-landscape-functions.jl")

Random.seed!(1234);
dims = (50,50);
S = 100;
W, P = fcn_spatial_weights(dims, 3);
X=10 .+rand(size(W,1),1);
Wp = fcn_spatial_AR(W, S; X = X, θ=0.9, σ=0.5);
x_hms = heatmap(reshape(X, dims));
hms = [heatmap(reshape(Wp[:,s], dims)) for s=1:4];
hms_mean = heatmap(reshape(mean(Wp,dims=2), dims));
plot(hms...;aspect_ratio = :equal, grid=false);
plot(hms_mean)