using HighDimMixedModels
using CSV
using DataFrames
using Lasso
using MLBase
using Random
using JLD2

# get data
data_dir = "data/real/gene_expressions"
ribo = CSV.read("$data_dir/riboflavingrouped.csv", DataFrame; 
    transpose = true, types = Dict(1 => String)) 
gene_expressions = Matrix{Float64}(ribo[:,3:end])
all_gene_names = names(ribo)[3:end]
y = ribo[:,2]
N = length(y)
# group info
grp = readlines("$data_dir/riboflavingrouped_structure.csv")
grp = [replace(x, "\"" => "") for x in grp]
# Control 
control = Control()
control.trace = true

# Fit random intercept model and save the names of the genes with non-zero effects in this model
# Note that the model results in no variation in intercepts 
# (so effectively the same as a LASSO without random effects)
X = ones(N, 1)
λ = 43
control.tol = 1e-3
Random.seed!(1234)
# we choose to standardize
int_fit = hdmm(X, gene_expressions, y, grp; 
    penalty="scad", λ=λ, ψstr="ident", control=control)
beta = int_fit.fixef
println("Number of non-zero coefs in initial fit: $(int_fit.nz), bic is $(int_fit.bic)")
idxs = findall(beta[2:end] .!= 0) #skip the intercept
gene_names = all_gene_names[idxs]
save_object("data/real/gene_expressions/gene_names.jld2", gene_names)


# cycle through each non-zero coefficient estimated from the random intercept fit (which resulted in an intercept variance of 0)
# for each such coefficient, associate a random effect to the corresponding predictor
# then fit the model with this random effect
res = Vector{Any}(undef, length(idxs))
for (i, idx) in enumerate(idxs)
    println("idx = $idx")
    # Include random effect for only the idx-th predictor (no random intercept!)
    Z = gene_expressions[:, [idx]]
    # standardize Z in line with Schelldorfer since algorithm itself doesn't standardize random effect design matrix
    Z = (Z .- mean(Z, dims=1)) ./ std(Z, dims=1)
    # Fixed effect design excludes the idx-th predictor
    G = gene_expressions[:, Not(idx)]
    # Unpenalized fixed effect design includes the intercept and the idx-th predictor
    local X = [ones(N) Z]
    # fit model
    local control = Control()
    control.trace = true
    control.tol = 1e-3
    local λ = 45
    Random.seed!(1234)
    est = hdmm(X, G, y, grp, Z; penalty="scad", λ=λ, ψstr="ident", control=control)
    res[i] = est
end
save_object("data/real/gene_expressions/gene_results.jld2", res)


###################################################################
# Fit a final model with random effects assigned to only the genes 
# with the highest random effects from the previous step
res = load_object("data/real/gene_expressions/gene_results.jld2")
gene_names = load_object("data/real/gene_expressions/gene_names.jld2")
# Extract ψs and find genes with high random effects
ψs = [x.ψ[1] for x in res]
kappa = 0.04
lre_mask = ψs .> kappa
println("Number of genes with high random effect: $(sum(lre_mask))")
println("Gene names with high random effect: $(gene_names[lre_mask])")
println("ψs for genes with high random effect: $(ψs[lre_mask])")
println("BICS for genes with high random effect: $(map(x -> x.bic, res[lre_mask]))")
rand_idx = indexin(gene_names[lre_mask], all_gene_names)
Z = gene_expressions[:, rand_idx]
# standardize Z in line with Schelldorfer since algorithm itself doesn't standardize random effect design matrix
Z = (Z .- mean(Z, dims=1)) ./ std(Z, dims=1)
X = ones(N, 1)
G = gene_expressions
λ = 43
control.tol = 1e-3
Random.seed!(1234)
@time final_fit = hdmm(X, G, y, grp, Z; penalty="scad", λ=λ, ψstr="diag", control=control)
save_object("data/real/gene_expressions/final_fit.jld2", final_fit)

# Inspection of final_fit
final_fit_loaded = load("data/real/gene_expressions/final_fit.jld2")["single_stored_object"]
print("Final fit: L is $(final_fit.L), σ² is $(final_fit.σ²)")
println("Number of non-zero coefs in final fit: $(final_fit.nz), bic is $(final_fit.bic)")
HighDimMixedModels.coeftable(final_fit, ["intercept"; all_gene_names])