using DataFrames

function ratios()

    ratio_matrix = DataFrame(μS=Float64[], μB=Float64[], νS=Float64[], νB=Float64[])

    for i in [1/2, 1, 2]
        for j in [1/2, 1, 2]
            for k in [1/2, 1, 2]
                for l in [1/2, 1, 2]    
                    push!(ratio_matrix, (μS=i, μB=j, νS=k, νB=l))
                end
            end
        end
    end

    n_rows = size(ratio_matrix, 1)

    constrained_matrix = DataFrame(μS=Float64[], μB=Float64[], νS=Float64[], νB=Float64[])

    for i in 1:n_rows

        μS = ratio_matrix[i, :μS]
        μB = ratio_matrix[i, :μB]
        νS = ratio_matrix[i, :νS]
        νB = ratio_matrix[i, :νB]

        if (
            (abs(μS/μB) < 3) && (abs(μS/νS) < 3) && (abs(νB/νS) < 3) &&
            (abs(μS/μB) > 1/3) && (abs(μS/νS) > 1/3) && (abs(νB/νS) > 1/3) &&
            (μS != 1 || μB != 1 || νS != 1 || νB != 1) && 
            (μB != 1/2 || νS != 1/2) && (μB != 1/2 || νB != 2)
        )
            push!(constrained_matrix, (μS=μS, μB=μB, νS=νS, νB=νB))
        end
    end

    return constrained_matrix
end

function resum_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order)
    ratio_matrix = ratios()
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        μS_ratio = ratio_matrix[i, :μS]
        μJ_ratio = ratio_matrix[i, :μB]
        νS_ratio = ratio_matrix[i, :νS]
        νJ_ratio = ratio_matrix[i, :νB]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                                μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, μS_ratio=μS_ratio, νS_ratio=νS_ratio)

        push!(delta_list, abs.(varied_list - central_list))
    end

    chi_length = length(χ_list)
    resum_list = []

    for i in 1:chi_length
        set = [delta_list[j][i] for j in 1:n_rows] 
        push!(resum_list, maximum(set)) 
    end

    return resum_list

end

function fixed_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order)
    ratio_matrix = DataFrame(μS=Float64[1/2,2], μB=Float64[1/2,2], νS=Float64[1/2,2], νB=Float64[1/2,2], μH=Float64[1/2,2])
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        μS_ratio = ratio_matrix[i, :μS]
        μJ_ratio = ratio_matrix[i, :μB]
        νS_ratio = ratio_matrix[i, :νS]
        νJ_ratio = ratio_matrix[i, :νB]
        μH_ratio = ratio_matrix[i, :μH]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                                μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, μS_ratio=μS_ratio, νS_ratio=νS_ratio, μH_ratio=μH_ratio)

        push!(delta_list, abs.(varied_list - central_list))
    end

    chi_length = length(χ_list)
    fixed_list = []

    for i in 1:chi_length
        set = [delta_list[j][i] for j in 1:n_rows] 
        push!(fixed_list, maximum(set)) 
    end

    return fixed_list

end

function NP_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order)
    ratio_matrix = DataFrame(bmax=[1/2,2])
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        bmax_ratio = ratio_matrix[i, :bmax]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                                bmax_ratio=bmax_ratio)

        push!(delta_list, abs.(varied_list - central_list))
    end

    chi_length = length(χ_list)
    NP_list = []

    for i in 1:chi_length
        set = [delta_list[j][i] for j in 1:n_rows] 
        push!(NP_list, maximum(set)) 
    end

    return NP_list

end