function CARA(w, α=0.0)
    u = α==0 ? w : (1 .-exp.(-α .* w)) ./ α;
    return u
end

function CRRA(w, θ=1.0)
    u = θ==1 ? log.(w) : w.^(1-θ) ./ (1-θ)
    return u
end

function EU(W, f, p = ones(size(W,1))/size(W,1))
    U = f.(W);
    EU = U' * p;
    return EU
end