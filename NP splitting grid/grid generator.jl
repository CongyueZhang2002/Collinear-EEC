using Distributed
addprocs(14)

@everywhere using CSV, DataFrames
@everywhere begin
    include("..\\core\\constants.jl")
    include("moments.jl")
    const Na1b = 3000
    const Na2  = 1000
    const a1bmin=0.0; const a1bmax=15.0; const p=2
    const a2min=0.9; const a2max=1.3
    const a1b_vals = a1bmin .+ (a1bmax - a1bmin) .* range(0, 1; length=Na1b).^p
    const a2_vals  = range(a2min, a2max; length=Na2)
    gamma_cols = [
        "γqq0","γqg0","γgq0","γgg0",
        "dγqq0","dγqg0","dγgq0","dγgg0",
        "ddγqq0","ddγqg0","ddγgq0","ddγgg0",
        "γqq1","γqg1","γgq1","γgg1",
        "dγqq1","dγqg1","dγgq1","dγgg1",
        "γqq2","γqg2","γgq2","γgg2"
    ]
    const γsyms = Symbol.(gamma_cols)
    const Ncols = length(γsyms)
    @inline linear_index(i_a1b, i_a2) = (i_a1b - 1) * Na2 + i_a2
end

using SharedArrays, ProgressMeter

const Ntot = Na1b * Na2
A1B = SharedArray{Float64}(Ntot)
A2  = SharedArray{Float64}(Ntot)
COL = [SharedArray{Float64}(Ntot) for _ in 1:Ncols]
threshold = 1e-8

prog_chan = RemoteChannel(() -> Channel{Nothing}(Na1b))
@everywhere prog_chan = $prog_chan

# Progress task (master)
@async begin
    prog = Progress(Na1b; desc="Building grid", showspeed=true)
    for _ in 1:Na1b
        take!(prog_chan)
        next!(prog)
        # Base.flush(stdout)  # Uncomment if VS Code still buffers
    end
end

@sync @distributed for i_a1b in 1:Na1b
    a1b = a1b_vals[i_a1b]
    for i_a2 in 1:Na2
        a2  = a2_vals[i_a2]
        idx = linear_index(i_a1b, i_a2)
        rawvals = buildgrid_func(a1b=a1b, a2=a2)
        A1B[idx] = a1b
        A2[idx]  = a2
        @inbounds for j in 1:Ncols
            v = rawvals[j]
            COL[j][idx] = (abs(v) < threshold) ? 0.0 : v
        end
    end
    put!(prog_chan, nothing)   # signal completion of one outer index
end

# Assemble / write
df = DataFrame(a1b = Array(A1B), a2 = Array(A2))
for (j,s) in enumerate(γsyms)
    df[!, s] = Array(COL[j])
end
outpath = raw"C:\Users\congyue zhang\Desktop\EEC collinear bspace new\NP splitting grid\splitting_moments_grid.csv"
CSV.write(outpath, df)
println("\nWrote CSV to: ", outpath)