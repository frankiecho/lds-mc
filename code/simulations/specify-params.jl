using CSV, Setfield

include("../../code/functions/type-defs.jl")

# Default parameters
param_default = LandscapeParameters((20,20), 100, 0.8, 6, 0.01, 20, 0.99)
param_vec = [param_default]
nsims = 5

# Vectors to test
dims_vec = []
budget_vec = [40, 60, 80, 100]
yy_vec = [300, 1000, 10000]
ρ_vec = [0, 0.2, 0.4, 0.6]
σ_vec = [2, 4, 8, 10]
η_vec = [0.1, 0.05, 0.001, 0.0001, 0.00001, 0] 
β_vec = [0.5, 0.7, 0.9, 0.95]

for budget in budget_vec
    local param_t = deepcopy(param_default)
    param_t.budget = budget
    push!(param_vec, param_t)
end

for yy in yy_vec
    local param_t = deepcopy(param_default)
    param_t.yy = yy
    push!(param_vec, param_t)
end

for ρ in ρ_vec
    local param_t = deepcopy(param_default)
    param_t.ρ = ρ
    push!(param_vec, param_t)
end

for σ in σ_vec
    local param_t = deepcopy(param_default)
    param_t.σ = σ
    push!(param_vec, param_t)
end

for η in η_vec
    local param_t = deepcopy(param_default)
    param_t.η = η
    push!(param_vec, param_t)
end

for dim in dims_vec
    local param_t = deepcopy(param_default)
    param_t.dims = dim
    push!(param_vec, param_t)
end

for β in β_vec
    local param_t = deepcopy(param_default)
    param_t.β = β
    push!(param_vec, param_t)
end

param_df = DataFrame(param_vec)
CSV.write("output/param_df_$(nsims).csv",param_df)