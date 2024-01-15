# Functions to simulate returns in a randomly generated landscape
using SparseArrays, LinearAlgebra, Random, Distributions, Revise, ThreadsX, DataFrames
include("optim-functions.jl")

function fcn_meshgrid(x,y)
    X = x' .* ones(y[end]);
    Y = ones(x[end])' .* y;
    return(X,Y)
end

function fcn_spatial_weights(dims = (5,5); boundary::Integer = 0, order = 1, queen = true, distance = false, bandwidth = 1, α=1)
    # Generates spatial weights matrices based on contiguity
    # dims: 2d vector of the number of x and y dimensions
    dims_boundary = dims .+ (boundary*2);
    x = 1:(dims_boundary[1]);
    y = 1:(dims_boundary[2]);
    X,Y = fcn_meshgrid(x,y);
    Xl = reshape(X, :, 1)';
    Yl = reshape(Y, :, 1)';
    P = spzeros(prod(dims_boundary), prod(dims_boundary));
    D = spzeros(prod(dims_boundary), prod(dims_boundary));
    for i = 1:prod(dims_boundary)
        if (distance) 
            dij = sqrt.((Xl.-Xl[i]).^2 .+ (Yl.-Yl[i]).^2);
            k = dij.^(-α); # Inverse distance
            k[i] = 0;
            k[dij .> bandwidth] .= 0;
            P[:,i] = k';
            D[:,i] = dij';
        elseif (queen)
            neighbour = (abs.(Xl.-Xl[i]) .≤ order) .* (abs.(Yl.-Yl[i]) .≤ order);
            P[:,i] = neighbour';
        else
            manhattan_dist = abs.(Xl.-Xl[i]) + abs.(Yl.-Yl[i]);
            P[:,i] = manhattan_dist' .≤ order;
        end
    end
    P[diagind(P)] .= 0;

    W = spzeros(size(P));
    for i = 1:prod(dims_boundary)
        W[i,:] = P[i,:] / sum(P[i,:]);
    end
    boundary_id = ((Xl .<= boundary) .| (Xl .> (dims[1]+boundary)) .| (Yl .<= boundary) .| (Yl .> (dims[2]+boundary))) |> vec;
    return W, P, boundary_id, D;
end

function fcn_spatial_AR(W, S=1; X=rand(size(W,1),1), β=1, ρ=0, γ=0, θ=0, σ=1)
    # Simulates S realisations of a spatially autoregressive data generating process
    #
    # SAR: (I-ρW)^-1 * (Xβ + WXγ) + (I-ρW)^-1 *(I-θW)^-1 ϵ
    # ϵ ~ N(0, σ^2)
    #
    # Arguments
    # W: row-normalised spatial weights matrix
    # S: number of realisations

    N = size(W,1);
    ϵ = Normal(0, σ^2);
    E = rand(ϵ,N,S);
    A = inv(Matrix(I - ρ*W)) * (X*β + W*X*γ);
    B = inv(Matrix(I - ρ*W))*inv(Matrix(I - θ*W));
    Y = repeat(A, 1, S) + B'*E;
    return Y;
end

function fcn_spatial_shock(W, S::Integer, n_shocks::Int, shock_size; p=ones(size(W,1))/size(W,1))
    N = size(W,1);
    V = zeros(N,S);
    NSS = zeros(S);
    sl = Array{Vector}(undef, S, 1)
    for s=1:S
        #nss = rand(n_shocks); # Sample number of shocks from distribution n_shocks
        #shock_loc = wsample(1:N, p, round(Int, nss); replace = false);

        # Draw using a uniform distribution of independent draws with underlying probability
        shock = rand(N) .< p;
        nss = sum(shock);
        shock_loc = findall(x -> x .> 0, shock)
        if nss > n_shocks
            shock_loc = sample(shock_loc, n_shocks; replace = false)
            nss = n_shocks
        end
        
        for i=shock_loc
            sss = round(Int, rand(shock_size)); # Sample shock size from distribution
            if (sss == 0)
                continue
            elseif (sss == 1)
                V[i,s] = 1;
                continue
            end

            w = vec(W[:,i]);
            subset_index = w .> 0;
            neighbour_index = (1:length(w))[subset_index];
            wv = w[subset_index];
            shock_cells = wsample(neighbour_index, wv, sss-1; replace = false);
            push!(shock_cells, i);
            V[shock_cells,s] .= 1;
        end
        NSS[s] = nss;
        sl[s] = shock_loc;
    end
    return V, NSS, sl;
end

function fcn_spatial_shock_gev(W, S::Integer; p=1/size(W,1), dist=Distributions.TDist(3), T::Integer=100)
    N = size(W,1);
    V = zeros(N,S);
    NSS = zeros(S);
    sl = Array{Vector}(undef, S, 1)
    ss = Array{Vector}(undef, S, 1)

    threshold = quantile(dist, 1-p); # Threshold of exceedance from inverse CDF of the specified distribution

    for s=1:S
        # Draw using a uniform distribution of independent draws with underlying probability
        ζ = rand(dist, T);
        nss = sum(ζ .> threshold);
        shock_loc = randperm(N)[1:nss]
        shock_size = Int.(ceil.((N/100)*(ζ[ζ .> threshold] .- threshold))); # Shock size relative to total landscape size
        
        for i=1:nss
            sss = min(shock_size[i], size(W,1)-1); # Sample shock size from distribution
            if (sss == 0)
                continue
            elseif (sss == 1)
                V[shock_loc[i],s] = 1;
                continue
            end

            w = vec(W[:,i]);
            subset_index = w .> 0;
            neighbour_index = (1:length(w))[subset_index];
            wv = w[subset_index];
            shock_cells = wsample(neighbour_index, wv, sss-1; replace = false);
            V[shock_cells,s] .= 1;
        end
        
        NSS[s] = nss;

    end
    return V, NSS;
end


function fcn_reshape_lds(D, dims)
    dims_round = round.(dims);
    return reshape(D, dims_round);
end

function fcn_generate_landscape(dims=(40,40); S=(prod(dims)+1), yy=0, ρ=rand()*0.99, σ=rand(Chisq(10)), η=1/prod(dims))
    # Generates a (n1 x n2) grid cells of landscape describing benefits of protecting each cell, with random shocks that wipe out returns in a given set of cells under every scenario.

    # dims: dimensions of the landscape grid
    # S: number of scenarios to generate
    # yy: quantity to add to average returns from the spatial AR process
    # n_shock: a statistical distribution describing the number of shocks in each scenarios
    # shock_size: a statistical distribution describing the number of cells affected by a shock

    #   Wps: benefits after shocks (for S scenarios)
    #   Wp: spatially-autocorrelated "base" benefits
    #   SS: shock locations
    #   pw_rs: probability of shock (spatially-autocorrelated)

    W, P, b = fcn_spatial_weights(dims; boundary = 1, distance = true, bandwidth = sqrt(sum(dims.^2)),  α=2);
    X = rand(Normal(0,1), size(W,1))
    Wp   = fcn_spatial_AR(W, S; X = X, ρ=ρ, σ= σ) .+ yy;
    Wp = Wp[b .== 0, :];
    Wp[isless.(Wp,0)] .= 0;

    if (sum(Wp .< 0) > 0)
        throw("Some elements in Wp are less than zero")
    end

    #pw   = fcn_spatial_AR(W, 1; X = rand(size(W,1),1), ρ=0.5, σ=5) |> vec; # Probability of shock
    #pw   = pw[b .== 0,:];
    #pw[pw .< 0] .= 0;
    #pw   = pw ./ sum(pw);
    #SS, nss, sl   = fcn_spatial_shock(W[b .== 0, b .== 0], S, n_shocks, shock_size, p=vec(pw));
    SS, nss   = fcn_spatial_shock_gev(W[b .== 0, b .== 0], S; p=η)
    Wps  = copy(Wp);
    Wps[isless.(1e-5,SS)] .= 0;
    #pw_rs = fcn_reshape_lds(pw, dims);

    return Landscape(dims, Wps, Wp, SS, W[b .== 0, b .== 0], nss);
end

function fcn_get_location_scale(R::Matrix, budget::Real=100, pct_range = (5,95)) 
    # Returns the location and scale parameters used for scaling returns for the utility functions
    ev_soln = fcn_optim_ev(R; budget = budget)
    ev_returns = R' * ev_soln
    location = percentile(ev_returns, pct_range[1])
    scale = percentile(ev_returns, pct_range[2]) - percentile(ev_returns, pct_range[1]);
    return location, scale
end

function fcn_get_w_random(R::Matrix, budget::Real=100, nsims::Real=100)
    # Gets returns of a set of random x in the set
    rand_soln = zeros(size(R, 1), nsims);

    for sim in range(1, nsims)
        s = randperm(size(R, 1))[1:budget];
        rand_soln[s, sim] .= 1;
    end

    return mean(R' * rand_soln)
end

function fcn_map_ef(R::Matrix, optim_func::Function, budget::Real=100, λ::AbstractVector=0:0.1:1, β::Real=0.9)
    # Maps a given function to identify the efficiency frontier across λ values
    solutions = map(l->optim_func(-R; budget = budget, λ=l, β=β),λ) |> e->mapreduce(permutedims, vcat, e) |> transpose
    RS = R' * solutions;
    ef = EfficiencyFrontier(solutions, RS, λ, optim_func);
    return ef
end

function fcn_evaluate_ef(ef::EfficiencyFrontier, u::UtilityFunction)
    # Evaluates solutions on the efficiency frontier based on the utility function
    return mapslices(u.CE, ef.returns; dims = 1)
end
