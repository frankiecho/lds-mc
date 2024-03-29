function CARA(w, α=0.0)
    u = α == 0 ? w : (1 .- exp.(-α .* w)) ./ α
    return u
end

function CARA_inv(u, α=0.0)
    if α == 0.0
        return u
    else
        return -log.(1 .- α .* u) ./ α
    end
end

function CRRA(w, θ=1.0)
    u = θ == 1.0 ? log.(w) : w .^ (1 - θ) ./ (1 - θ)
    return u
end

function CRRA_norm(w, θ=1.0, a=2.0, b=1.0)
    return (CRRA(w, θ) .- CRRA(b, θ)) ./ (CRRA(a, θ) .- CRRA(b, θ))
end

function CRRA_inv(u, θ=1.0)
    w = θ == 1.0 ? exp.(u) : ((1 - θ) .* u) .^ (1 ./ (1 - θ))
    return w
end

function CRRA_inv_norm(u, θ=1.0, a=2.0, b=1.0)
    CRRA_inv((CRRA(a, θ) - CRRA(b, θ)) * u .+ CRRA(b, θ), θ)
end

function piecewise_linear(c, G, r, q)
    # c: consumption
    # G: goal (breakpoint of the piecewise linear utility function)
    # r: cost for not achieving goal by each unit
    # q: value of additional gains over G
    function fcn_piecewise_elementwise(c, G, r, q)
        if c > G
            return q * (c - G) + G
        else
            return r * (c - G) + G
        end
    end
    return fcn_piecewise_elementwise.(c, G, r, q)
end

function EU(W, f, p=ones(size(W, 1)) / size(W, 1))
    U = f.(W)
    EU = U' * p
    return EU
end

