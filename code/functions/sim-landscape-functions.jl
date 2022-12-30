# Functions to simulate returns in a randomly generated landscape
using SparseArrays, LinearAlgebra, Plots, Random, Distributions
function fcn_meshgrid(x,y)
    X = x' .* ones(y[end]);
    Y = ones(x[end])' .* y;
    return(X,Y)
end

function fcn_spatial_weights(dims = (5,5), order = 1, queen = true)
    # Generates spatial weights matrices based on contiguity
    # dims: 2d vector of the number of x and y dimensions

    x = 1:dims[1];
    y = 1:dims[2];
    X,Y = fcn_meshgrid(x,y);
    Xl = reshape(X, :, 1)';
    Yl = reshape(Y, :, 1)';
    P = spzeros(prod(dims), prod(dims));

    for i = 1:prod(dims)
        if (queen)
            neighbour = (abs.(Xl.-Xl[i]) .≤ order) .* (abs.(Yl.-Yl[i]) .≤ order);
            P[:,i] = neighbour';
        else
            manhattan_dist = abs.(Xl.-Xl[i]) + abs.(Yl.-Yl[i]);
            P[:,i] = manhattan_dist' .≤ order;
        end
    end
    P[diagind(P)] .= 0;

    W = spzeros(size(P));
    for i = 1:prod(dims)
        W[i,:] = P[i,:] / sum(P[i,:]);
    end
    boundary = isequal.(Xl, 1) .| isequal.(Xl, dims[1]) .| isequal.(Yl, 1) .| isequal.(Yl, dims[2])
    return W, P, boundary
end

function fcn_spatial_AR(W, S=1; X=rand(size(W,1),1), β=1, ρ=0, γ=0, θ=0, σ=1)
    N = size(W,1);
    ϵ = Normal(0, σ);
    E = rand(ϵ,N,S);
    Y = zeros(N,S);
    A = inv(Matrix(I - ρ*W)) * (X*β + W*X*γ);
    B = inv(Matrix(I - ρ*W))*inv(Matrix(I - θ*W));
    Y = repeat(A, 1, S) +B'*E;
    return Y;
end

