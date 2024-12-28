using ZipFile
using CSV
using DataFrames
using HighDimMixedModels
using Serialization


function extract_gwas_data(data)
    X_names = [col for col in names(data) if startswith(col, "X")]
    G_names = [col for col in names(data) if startswith(col, "G")]
    grp = string.(data[:, 1])
    X = Matrix{Float64}(data[:, X_names])
    G = Matrix{Float64}(data[:, G_names])
    y = data[:, end]
    return X, G, y, grp 
end

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

function slim_model(est::HDMModel)
    return SlimHDMModel(
        est.init_coef.βstart,
        est.log_like,
        est.aic,
        est.bic,
        est.nz,
        est.fixef,
        est.ψ,
        est.σ²
        )
end

function run_model(data, lambdas, penalty, cov_structure, tol)
    X, G, y, grp = extract_gwas_data(data)
    control = Control()
    control.trace = true
    control.tol = tol
    results = []
    # Fit model over range of lambdas
    for λ in lambdas
        println("λ is $λ")
        est = hdmm(X, G, y, grp;
            penalty=penalty, λ=λ, ψstr=cov_structure, control=control)
        if est.nz < 50
            display(est)
        end
        slim_est = slim_model(est)
        push!(results, (λ, slim_est))
    end
    return results
end

function process_experiment(zip_path, penalty, cov_structure, lambdas, tol)
    zip_path = "../../data/GWAS/$(zip_path)"
    reader = ZipFile.Reader(zip_path)
    results = []
    for (i, file) in enumerate(reader.files)
        if endswith(file.name, ".csv")
            println("Processing file $i: $(Base.basename(file.name))")
            data = CSV.read(file, DataFrame)
            # Run model for each lambda 
            model_results = run_model(data, lambdas, penalty, cov_structure, tol)
            for (λ, slim_est) in model_results
                push!(results, (Base.basename(file.name), λ, slim_est))
            end
        end
    end
    close(reader)

    # remove .zip extension
    setting = split(Base.basename(zip_path), ".")[1]
    serialize("$(setting)_$(penalty)-newresults.txt", results)
end

function main()
    open("settings.txt", "r") do file
        for line in eachline(file)
            parts = split(strip(line), ", ")
            zip_path = String(parts[1])
            penalty = String(parts[2])
            cov_structure = String(parts[3])
            min_lambda = parse(Int, parts[4])
            max_lambda = parse(Int, parts[5])
            incr_lambda = parse(Int, parts[6])
            tol = parse(Float64, parts[7])
            lambdas = range(min_lambda, step=incr_lambda, stop=max_lambda)
            process_experiment(zip_path, penalty, cov_structure, lambdas, tol)
        end
    end
end

main()
