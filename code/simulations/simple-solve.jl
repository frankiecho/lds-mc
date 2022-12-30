using DataFrames
using Distributions
using JuMP
using LinearAlgebra
using Random
using Gurobi
using Revise
using StatsPlots
using Statistics
using StatsBase

includet("../functions/optim-functions.jl");

# Cost distribution
n = 100;
K = 10000;
mu = -rand(n);
sigma = Diagonal(rand(n));
dist = MvNormal(mu, sigma);
Y = rand(dist, K);

β = 0.9;
p = ones(K);
p = p/sum(p);
budget = 3;

ev_soln = fcn_optim_cvar(Y; p = p, budget = 10, β = 0.9, λ = 0);
mix_soln = fcn_optim_cvar(Y; p = p, budget = 10, β = 0.9, λ = 0.6);
cvar_soln = fcn_optim_cvar(Y; p = p, budget = 10, β = 0.9, λ = 1);
mv_soln = fcn_optim_mv(Y; budget = 10, λ = 1);

density(Y' * ev_soln)
density!(Y' * cvar_soln)
density!(Y' * mix_soln)
density!(Y' * mv_soln)

w = ProbabilityWeights(p);



