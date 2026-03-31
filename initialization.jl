using Distributed

@everywhere begin
    using CSV, DataFrames, Interpolations, StaticArrays

    const Sigma_names = [
        "ΣbLL9","ΣbNLL9","ΣbNNLL9",
        "ΣbLL8","ΣbNLL8","ΣbNNLL8",
        "ΣbLL7","ΣbNLL7","ΣbNNLL7",
        "ΣbLL6","ΣbNLL6","ΣbNNLL6",
        "ΣbLL5","ΣbNLL5","ΣbNNLL5",
        "ΣbLL4","ΣbNLL4","ΣbNNLL4",
        "ΣbLL3","ΣbNLL3","ΣbNNLL3",
        "ΣbLL2","ΣbNLL2","ΣbNNLL2",
        "ΣbLL1","ΣbNLL1","ΣbNNLL1",
        "ΣbLL0","ΣbNLL0","ΣbNNLL0"
    ]
    const NΣ = length(Sigma_names)

    const a1b_vals = sort!(unique(df_Sigma_coeffs.a1b))
    const a2_vals  = sort!(unique(df_Sigma_coeffs.a2))
    const Na1b = length(a1b_vals);  const Na2 = length(a2_vals)
    @assert size(df_Sigma_coeffs,1) == Na1b * Na2

    # Build matrices (Na1b, Na2)
    mats = Matrix{Float64}[]
    for name in Sigma_names
        push!(mats, reshape(df_Sigma_coeffs[!, name], Na2, Na1b)' )
    end

    # Pack into Array of SVector{NΣ}
    grid = Array{SVector{NΣ,Float64}}(undef, Na1b, Na2)
    @inbounds for i in 1:Na1b, j in 1:Na2
        grid[i,j] = SVector{NΣ,Float64}( m[i,j] for m in mats )
    end

    const Sigma_itp = interpolate((a1b_vals, a2_vals), grid, Gridded(Linear()))

    # Free DataFrame if not needed further
    df_Sigma_coeffs = nothing

    # ---- Public API ----

    # Positional + keyword variants returning SVector
    @inline Sigma_coeffs_func(a1b::Float64, a2::Float64) = Sigma_itp(a1b, a2)
    @inline Sigma_coeffs_func(; a1b::Float64, a2::Float64) = Sigma_itp(a1b, a2)

    # Allocate a fresh Vector (slowest of the three, but matches old type)
    #@inline Sigma_coeffs_vec(a1b::Float64, a2::Float64) = collect(Sigma_itp(a1b, a2))
    #@inline Sigma_coeffs_vec(; a1b::Float64, a2::Float64) = collect(Sigma_itp(a1b, a2))

    # In-place fill into an existing Vector{Float64}
    #@inline function Sigma_coeffs_func!(buf::Vector{Float64}, a1b::Float64, a2::Float64)
    #    @assert length(buf) == NΣ
    #    v = Sigma_itp(a1b, a2)
    #    @inbounds @simd for k in 1:NΣ
    #        buf[k] = v[k]
    #    end
    #    return buf
    #end
    #@inline Sigma_coeffs_func!(buf::Vector{Float64}; a1b::Float64, a2::Float64) =
    #    Sigma_coeffs_func!(buf, a1b, a2)
end
