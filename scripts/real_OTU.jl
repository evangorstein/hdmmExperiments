using HighDimMixedModels
using StatsBase
using JLD2
using CSV
using DataFrames
using Plots
using RCall

# get data
# note that response y is square root of age and the design matrix is library size normalized
data_dir = "data/real/OTU"
d = CSV.read("$data_dir/normalized_data.csv", DataFrame)
otu_taxonomies = CSV.read("$data_dir/final_taxonomies.csv", DataFrame)
N = size(d)[1]
grp = convert(Vector{String}, d.country)
# USA: 1, Malawi: 2, Venezuela: 3
countmap(grp) 
otus = Matrix(d[:,3:(end-1)])
y = d.sqrt_age
X = ones(N, 1)

control = Control()
control.trace = 3
log_y = log.(y)
int_fit_orig = hdmm(
X, otus, log_y, grp;
penalty="scad", λ=200, ψstr="ident", control=control
)

# check diagnostics
scatter(int_fit_orig.fitted, int_fit_orig.resid)
scatter(otus[:,365], log_y)
mask_outly = int_fit_orig.resid .> 1.5
scatter(int_fit_orig.fitted, log_y)
annotate!(
int_fit_orig.fitted[mask_outly].+0.01, 
log.(y)[mask_outly], 
text.(d.sample_id[mask_outly], :red, :left, 11)
)


# fit model with dropped otus and with quadratic for 365th otu
otus = [otus otus[:,365].^2] 
@time int_fit = hdmm(
X[.!mask_outly,:], otus[.!mask_outly,:], log_y[.!mask_outly], grp[.!mask_outly]; 
penalty="scad", λ=200, ψstr="ident", control=control
)
scatter(int_fit.fitted, int_fit.resid)
scatter(int_fit.fitted, log_y[.!mask_outly])
save_object("data/real/OTU/otu_fit.jld2", int_fit)



# load otu_fit.jld2
otu_fit = load_object("data/real/OTU/otu_fit.jld2")



# get indices of non-zero coefficients besides intercept in the initial estimates
init_effects = otu_fit.init_coef.βstart[2:end-1]
idx_nz_init = findall(!iszero, init_effects)
# get taxonomies of non-zero coefficients and their values
coef_nz_init = [otu_taxonomies[idx_nz_init,end-3:end] init_effects[idx_nz_init]] 
# get max and min absolute value of non-zero coefficients
max_coef = maximum(abs.(coef_nz_init[:,end]))
min_coef = minimum(abs.(coef_nz_init[:,end]))

# get indices of non-zero coefficients besides intercept in the final estimates
effects = otu_fit.fixef[2:end-1]
idx_nz = findall(!iszero, effects)
# get names of non-zero coefficients and their values
coef_nz = [idx_nz otu_taxonomies[idx_nz, (end-2):end] ]
coef_nz.effect_estimates = effects[idx_nz]
rename!(coef_nz, "x1" => "OTU_index")
R"
library(tinytable)
$coef_nz |>
    tt(digits = 2) |>
    style_tt(align = \"c\") |>
    print(\"latex\")
"

#Calculate variance explained
raw_var = var(log_y)
resid_var = var(otu_fit.resid)
# same as otu_fit.σ²
var_exp = 1 - resid_var/raw_var
