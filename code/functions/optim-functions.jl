
using Gurobi, JuMP, Statistics, StatsBase, Distributed

# Optimisation of conditional value-at-risk with JuMP
function fcn_optim_cvar(Y; p = ones(size(Y,2),1)/size(Y,2), budget = 1, β = 0.9, λ::Float64 = 1.0)
    # Optimizes the conditional value-at-risk of a portfolio with N assets
    # Y: a N by K matrix, with K realisations of the uncertain vector of costs
    # p: probability of K events
    # budget: assets must sum to budget
    # β: CVaR tail quantile (probability that costs are lower than value-at-risk)
    # λ: risk-aversion parameter for the objective function with EV and CVaR optimisation

    N, K = size(Y);
    if ~isnothing(p)
        p = ones(K,1)/K;
    end
    μ = Y*p;

    model = Model(Gurobi.Optimizer);
    set_silent(model);
    @variable(model,0 <= x[1:N] <= 1);
    @variable(model, α);
    @variable(model, u[1:K] >= 0);
    @objective(model, Min, (1-λ)*sum(μ' * x) + λ*(α + sum(p'*u)/(1-β)))
    for k=1:K
        @constraint(model, -x' * Y[:,k] + α + u[k] >= 0);
    end
    @constraint(model, sum(x) == budget);
    optimize!(model);
    return(value.(x));
end

function fcn_optim_mv(Y; budget = 1, λ = 1)
    # Optimisation of mean-variance portfolio
    N, K = size(Y);
    Q = StatsBase.cov(Y'); # Weighted covariance
    μ = mean(Y);
    portfolio = Model(Gurobi.Optimizer);
    set_silent(portfolio)
    @variable(portfolio, x[1:N] >= 0)
    @objective(portfolio, Min, (1-λ)*sum(μ' * x) + λ * x' * Q * x)
    @constraint(portfolio, sum(x) == budget)
    optimize!(portfolio)
    return(value.(x));
end

# Optimisation of EV with JuMP
function fcn_optim_ev(Y; p = ones(size(Y,2),1)/size(Y,2), budget = 1)
    # Optimizes the conditional value-at-risk of a portfolio with N assets
    # Y: a N by K matrix, with K realisations of the uncertain vector of costs
    # p: probability of K events
    # budget: assets must sum to budget
    # β: CVaR tail quantile (probability that costs are lower than value-at-risk)
    # λ: risk-aversion parameter for the objective function with EV and CVaR optimisation

    N, K = size(Y);
    if ~isnothing(p)
        p = ones(K,1)/K;
    end
    μ = Y*p;

    model = Model(Gurobi.Optimizer);
    set_silent(model);
    @variable(model,0 <= x[1:N] <= 1);
    @objective(model, Min, sum(μ' * x))
    @constraint(model, sum(x) == budget);
    optimize!(model);
    return(value.(x));
end