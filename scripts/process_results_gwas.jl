using HighDimMixedModels
using Serialization
using CSV
using DataFrames

## The results for the GWAS experiments are stored in a slightly different format
## from those of the gene expression and OTUs. 

struct SlimHDMModel
    init_coef::Vector{Float64}
    log_like::Float64
    aic::Float64
    bic::Float64
    nz::Float64
    fixef::Vector{Float64}
    ψ::Matrix{Float64}
    σ²::Float64
end


data_files = ["data" * string(i) * ".csv" for i in 1:100]
coef_cols =  9 + 25 + 1001 # 9 model fit stats, 25 psi parameters, 1001 beta parameters
n_rows = 2 * 10 * 100 # top two results for each setting, data-set pair
best_results = Matrix{Any}(undef, n_rows, coef_cols)
i = 1 #Index of next row to fill in in best_results

dir_path = "sim_results/gwas_new/"
setting_names = readdir(dir_path)
setting_names = setting_names[endswith.(setting_names, "results.txt")]
# How many true non-zero parameters are there? They all come at the beginning of the vector
true_nz = 10 

#Iterate through each setting
for setting_name in setting_names
    
    # get serialized results for this setting
    setting = open(dir_path*setting_name, "r")
    setting_results = deserialize(setting)
    close(setting)

    #Iterate through each data file fit under this setting
    for file in data_files
        
        data_id = replace(file, ".csv" => "")
        data_results = filter(x -> x[1] == file, setting_results)
        
        # Get results from best two λs
        data_sorted_results = sort(data_results, by = x -> x[3].bic)
        data_best_results = data_sorted_results[1:2]
        
        #Add a single row to the best_results matrix for each row in data_best_results
        for result in data_best_results

            # Get number of true positives
            nz_ind = findall(result[3].fixef .!= 0)
            tp = length(intersect(nz_ind, 1:true_nz))

            # Remove suffix from setting_name
            setting_name = replace(setting_name, ".txt" => "")

            # Create row to get added to matrix of best results
            added_row = [setting_name, data_id, 
                        result[2], result[3].log_like, 
                        result[3].aic, result[3].bic, 
                        result[3].nz, tp, sqrt(result[3].σ²)]
            
            # Add psi
            vec_ψ = vec(result[3].ψ)
            added_row = [added_row; vec_ψ; fill(missing, 25-length(vec_ψ))]

            # Add beta estimates
            added_row = [added_row; result[3].fixef]

            # Add row to matrix of best results
            best_results[i,:] = added_row
            global i += 1
        end
    end
end

#Create dataframe
colnames = ["setting", "data_id", "lambda", "loglike", "aic", "bic", "n_nz", "tp", "sigma"]
colnames = vcat(colnames, ["psi_$i" for i in 1:25]...)
colnames = vcat(colnames, ["beta_$i" for i in 1:1001]...)

best_results_df = DataFrame(best_results[1:i-1,:], colnames)

# Define the output file name 
file_path = "sim_results/gwas_new/best_results.csv" 

# Write the data to file
CSV.write(file_path, best_results_df)
