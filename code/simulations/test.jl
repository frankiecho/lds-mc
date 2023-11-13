using Distributions, Plots
S = 160000
N = 1600;
p = 0.0001;
dist=Distributions.TDist(3);
threshold = quantile(dist, 1-p);
ζ = rand(dist, S);
shock_states = findall(x -> x .> threshold, ζ)
shock_size = Int.(ceil.((N/100).*(ζ[shock_states].-threshold)));

histogram(shock_size)