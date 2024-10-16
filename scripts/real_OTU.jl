using HighDimMixedModels
using RCall
using CSV
using DataFrames
using Lasso
using MLBase
using Random
using JLD2

# get data
# note that response y is square root of age and the design matrix is library size normalized
data_dir = "data/real/OTU"
d = CSV.read("$data_dir/tss_normalized_data.csv", DataFrame)
otu_taxonomies = CSV.read("$data_dir/taxonomies.csv", DataFrame)
N = size(d)[1]
grp = d.country
# USA: 1, Malawi: 2, Venezuela: 3
countmap(grp) 
otus = Matrix(d[:,2:(end-1)])
y = d.age
X = ones(N, 1)

control = Control()
control.trace = 3
λ = 75
control.tol = 1e-3
Random.seed!(1234)
# We choose to standardize because the data has been library sum scaled and so is very small numbers
@time int_fit = hdmm(X, otus, y, grp; 
    penalty="scad", λ=λ, ψstr="ident", control=control)

save_object("data/real/OTU/otu_fit.jld2", int_fit)


# load otu_fit.jld2
otu_fit = load_object("data/real/OTU/otu_fit.jld2")



# get indices of non-zero coefficients besides intercept in the initial estimates
init_effects = otu_fit.init_coef.βstart[2:end]
idx_nz_init = findall(!iszero, init_effects)
# get taxonomies of non-zero coefficients and their values
coef_nz_init = [otu_taxonomies[idx_nz_init,end-3:end] init_effects[idx_nz_init]] 
# get max and min absolute value of non-zero coefficients
max_coef = maximum(abs.(coef_nz_init[:,end]))
min_coef = minimum(abs.(coef_nz_init[:,end]))

# get indices of non-zero coefficients besides intercept in the final estimates
effects = otu_fit.fixef[2:end]
idx_nz = findall(!iszero, effects)
# get names of non-zero coefficients and their values
coef_nz = [idx_nz otu_taxonomies[idx_nz, (end-2):end] ]
coef_nz.effect_estimates = effects[idx_nz]
rename!(coef_nz, "x1" => "OTU_index")
R"
$coef_nz |>
    tt(digits = 2) |>
    style_tt(align = \"c\") |>
    print(\"latex\")
"

#Calculate variance explained
raw_var = var(y)
resid_var = var(otu_fit.resid)
# same as otu_fit.σ²
var_exp = 1 - resid_var/raw_var
