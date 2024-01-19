## Stores type definitions

using Revise, DataFrames
include("expected-utility-functions.jl")

struct Landscape
    # Two dimensional landscape over N scenarios
    dims::Tuple # Dimensions (n1 x n2)
    R::Matrix # Returns (including shock)
    RV::Matrix # Returns (without shock)
    SS::Matrix # Shock location (binary)
    W::Matrix # Spatial weights
    nss::Vector # Number of shocks
end

struct UtilityFunction
    # Defines a von-Neumann Morgenstern utility function that scales inputs based on location and scale
    type::String
    α::Real
    location::Real
    scale::Real
    U::Function
    inv::Function
    EU::Function
    CE::Function

    function UtilityFunction(type, α, location=1.0, scale=1.0)
        fcn_scale = w -> (w .- location) ./ scale
        fcn_scale_inv = w -> (w .* scale) .+ location
        if cmp(type, "CARA") == 0
            U = (w) -> (w |> fcn_scale |> w -> CARA(w, α))
            inv = (w) -> CARA_inv(w, α) |> fcn_scale_inv
        elseif cmp(type, "CRRA") == 0
            U = (w) -> CRRA(fcn_scale(w), α)
            inv = (w) -> fcn_scale_inv(CRRA_inv(w, α))
        else
            U = (w) -> piecewise_linear(w |> fcn_scale, 0.0, α, 1.0)
            inv = (w) -> w |> fcn_scale_inv
        end
        eu = (w) -> EU(w, U)
        ce = (w) -> inv(EU(w, U))
        return new(type, α, location, scale, U, inv, eu, ce)
    end
end

struct EfficiencyFrontier
    solutions::AbstractArray
    returns::AbstractArray # Returns for each solution under each scenario
    λ::AbstractVector
    f::Function
end

struct Result
    ev::Any
    ev_rv::Any # EV optimisation without shocks
    cvar::Any
    mstd::Any
end

struct MCResult
    ef::Result # Efficiency Frontier
    ce::Result # Certainty equivalents
    ce_max::Result # Certainty equivalents (with maximised λ)
    L::Landscape # Landscape
end

mutable struct LandscapeParameters
    dims::Tuple # Dimensions (n1 x n2)
    yy::Real # Base return
    ρ::Real # Spatial correlation
    σ::Real # Standard deviation of SAR process
    η::Real # risk threshold
    budget::Real # Budget
    β::Real
    shocks::Tuple
end