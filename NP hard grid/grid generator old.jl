using CSV, DataFrames, ProgressMeter

include("..\\core\\constants.jl")
include("hard.jl")
include("harmonic polylogs.jl")


# ── 1) Grid parameters
Na1b, Na2       = 2500, 400
a1bmin, a1bmax = 0.0, 25.0
a2min, a2max   = 0.5, 2.0
p               = 2

# ── 2) Generate spaced vectors
# Power-law spacing for a1b
a1b_vals = a1bmin .+ (a1bmax - a1bmin) .* collect(range(0, 1; length=Na1b)).^p
# Linear spacing for a2
a2_vals  = collect(range(a2min, a2max; length=Na2))

# ── 3) Column definitions
hard_cols = ["hq0","hq1","hq2",
            "hg0","hg1","hg2"]

# ── 4) Initialize empty DataFrame
col_syms = Symbol.(vcat(["a1b", "a2"], hard_cols))
df = DataFrame()
for sym in col_syms
    df[!, sym] = Float64[]
end

# ── 5) Fill DataFrame with progress meter and threshold small values
@showprogress for a1b in a1b_vals, a2 in a2_vals
    # call your Julia function
    rawvals = hard_buildgrid_func(a1b=a1b, a2=a2)
    # convert to array and zero out tiny entries
    vals = [abs(x) < 1e-8 ? 0.0 : x for x in rawvals]

    # assemble row and push
    row  = (; a1b=a1b, a2=a2, Pair.(Symbol.(hard_cols), vals)...)  
    push!(df, row)
end

# ── 6) Write out to CSV
outpath = raw"C:\Users\congyue zhang\Desktop\EEC collinear bspace new\NP hard grid\hard_moments_grid.csv"
CSV.write(outpath, df)
println("Wrote CSV to: ", outpath)