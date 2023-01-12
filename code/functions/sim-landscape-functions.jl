# Functions to simulate returns in a randomly generated landscape
using SparseArrays, LinearAlgebra, Plots, Random, Distributions
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
    for i = 1:prod(dims_boundary)
        if (distance) 
            dij = sqrt.((Xl.-Xl[i]).^2 .+ (Yl.-Yl[i]).^2);
            k = dij.^(-α); # Inverse distance
            k[i] = 0;
            k[dij .> bandwidth] .= 0;
            P[:,i] = k';
            
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
    return W, P, boundary_id;
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

function fcn_spatial_shock(W, S::Integer, n_shocks, shock_size; p=ones(size(W,1))/size(W,1))
    N = size(W,1);
    V = zeros(N,S);
    NSS = zeros(S);
    sl = Array{Vector}(undef, S, 1)
    for s=1:S
        nss = rand(n_shocks); # Sample number of shocks from distribution n_shocks
        shock_loc = wsample(1:N, p, nss; replace = false);
        for i=shock_loc
            sss = rand(shock_size); # Sample shock size from distribution
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

function fcn_reshape_lds(D, dims)
    dims_round = round.(dims);
    return reshape(D, dims_round);
end

function fcn_generate_landscape(dims, S, n_shocks, shock_size)