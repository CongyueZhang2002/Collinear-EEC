using CSV, DataFrames, ProgressMeter
        
df_splitting_moments = CSV.read("C:\\Users\\congyue zhang\\Desktop\\EEC collinear bspace new\\NP splitting grid\\splitting_moments_grid.csv", DataFrame)
df_hard_moments = CSV.read("C:\\Users\\congyue zhang\\Desktop\\EEC collinear bspace new\\NP hard grid\\hard_moments_grid.csv", DataFrame)     

include("Sigma.jl")
include("initialization.jl")

# ── 0) Grid parameters

const μJ_ratio = 1/1.9

# ── 1) Grid parameters
Na1b, Na2       = 1000, 1000
a1bmin, a1bmax  = 0.0, 15.0
a2min, a2max    = 0.9, 1.3
p               = 2

# ── 2) Generate spaced vectors
# Power-law spacing for a1b
a1b_vals = a1bmin .+ (a1bmax - a1bmin) .* collect(range(0, 1; length=Na1b)).^p
# Linear spacing for a2
a2_vals  = collect(range(a2min, a2max; length=Na2))

# ── 3) Column definitions
Sigma_cols = [
    "ΣbLL9",  "ΣbNLL9",  "ΣbNNLL9",
    "ΣbLL8",  "ΣbNLL8",  "ΣbNNLL8",
    "ΣbLL7",  "ΣbNLL7",  "ΣbNNLL7",
    "ΣbLL6",  "ΣbNLL6",  "ΣbNNLL6",
    "ΣbLL5",  "ΣbNLL5",  "ΣbNNLL5",
    "ΣbLL4",  "ΣbNLL4",  "ΣbNNLL4",
    "ΣbLL3",  "ΣbNLL3",  "ΣbNNLL3",
    "ΣbLL2",  "ΣbNLL2",  "ΣbNNLL2",
    "ΣbLL1",  "ΣbNLL1",  "ΣbNNLL1",
    "ΣbLL0",  "ΣbNLL0",  "ΣbNNLL0"
]

# ── 4) Initialize empty DataFrame
col_syms = Symbol.(vcat(["a1b", "a2"], Sigma_cols))
df = DataFrame()
for sym in col_syms
    df[!, sym] = Float64[]
end

# ── 5) Fill DataFrame with progress meter and threshold small values
@showprogress for a1b in a1b_vals, a2 in a2_vals
    # call your Julia function
    rawvals = Sigma_buildgrid_func(a1b=a1b, a2=a2)
    # convert to array and zero out tiny entries
    vals = [abs(x) < 1e-8 ? 0.0 : x for x in rawvals]

    # assemble row and push
    row  = (; a1b=a1b, a2=a2, Pair.(Symbol.(Sigma_cols), vals)...)  
    push!(df, row)
end

for nm in names(df)
    df[!, nm] = Float32.(df[!, nm])
end

# ── 6) Write out to CSV
outpath = raw"C:\Users\congyue zhang\Desktop\EEC collinear bspace new\NP Sigma grid error\Sigma_coeffs_grid_1d19.csv"
CSV.write(outpath, df)
println("Wrote CSV to: ", outpath)