using CSV, DataFrames, ProgressMeter

include("..\\core\\constants.jl")
include("moments.jl")

# ── 1) Grid parameters
Na1b, Na2       = 3000, 1000
a1bmin, a1bmax  = 0.0, 15.0
a2min, a2max    = 0.9, 1.3
p               = 2

# ── 2) Generate spaced vectors
# Power-law spacing for a1b
a1b_vals = a1bmin .+ (a1bmax - a1bmin) .* collect(range(0, 1; length=Na1b)).^p
# Linear spacing for a2
a2_vals  = collect(range(a2min, a2max; length=Na2))

# ── 3) Column definitions
gamma_cols = ["γqq0","γqg0","γgq0","γgg0",
              "dγqq0","dγqg0","dγgq0","dγgg0",
              "ddγqq0","ddγqg0","ddγgq0","ddγgg0",
              "γqq1","γqg1","γgq1","γgg1",
              "dγqq1","dγqg1","dγgq1","dγgg1",
              "γqq2","γqg2","γqgq2","γgg2"]

# ── 4) Initialize empty DataFrame
col_syms = Symbol.(vcat(["a1b", "a2"], gamma_cols))
df = DataFrame()
for sym in col_syms
    df[!, sym] = Float64[]
end

# ── 5) Fill DataFrame with progress meter and threshold small values
@showprogress for a1b in a1b_vals, a2 in a2_vals
    # call your Julia function
    rawvals = buildgrid_func(a1b=a1b, a2=a2)
    # convert to array and zero out tiny entries
    vals = [abs(x) < 1e-8 ? 0.0 : x for x in rawvals]

    # assemble row and push
    row  = (; a1b=a1b, a2=a2, Pair.(Symbol.(gamma_cols), vals)...)  
    push!(df, row)
end

# ── 6) Write out to CSV
outpath = raw"C:\Users\congyue zhang\Desktop\EEC collinear bspace broken\NP splitting grid\moments_grid.csv"
CSV.write(outpath, df)
println("Wrote CSV to: ", outpath)