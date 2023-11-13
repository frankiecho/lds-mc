using CSV, Setfield

include("../../code/functions/type-defs.jl")

# Default parameters
param_default = LandscapeParameters((40,40), 0, 0.8, 6, 0.0001)
param_vec = [param_default]

# Vectors to test
yy_vec = []
ρ_vec = [0, 0.2, 0.4, 0.6]
σ_vec = [0, 2, 4, 8, 10]
η_vec = [0.01, 0.001, 0.0001, 0.00001, 0]

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


param_df = DataFrame(param_vec)
CSV.write("output/param_df.csv",param_df)