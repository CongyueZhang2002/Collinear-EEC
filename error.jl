# https://arxiv.org/pdf/1307.1808

using DataFrames
using Distributions

# Band Method
function ratios()

    ratio_matrix = DataFrame(μS_ratio=Float64[], μJ_ratio=Float64[], νS_ratio=Float64[], νJ_ratio=Float64[])

    for i in [1/2, 1, 2]
        for j in [1/2, 1, 2]
            for k in [1/2, 1, 2]
                for l in [1/2, 1, 2]    
                    push!(ratio_matrix, (μS_ratio=i, μJ_ratio=j, νS_ratio=k, νJ_ratio=l))
                end
            end
        end
    end

    n_rows = size(ratio_matrix, 1)

    constrained_matrix = DataFrame(μS_ratio=Float64[], μJ_ratio=Float64[], νS_ratio=Float64[], νJ_ratio=Float64[])

    for i in 1:n_rows

        μS_ratio = ratio_matrix[i, :μS_ratio]
        μJ_ratio = ratio_matrix[i, :μJ_ratio]
        νS_ratio = ratio_matrix[i, :νS_ratio]
        νJ_ratio = ratio_matrix[i, :νJ_ratio]

        if (
            (νJ_ratio/νS_ratio <= 2) && (νJ_ratio/νS_ratio >= 1/2) &&
            (μJ_ratio/μS_ratio <= 2) && (μJ_ratio/μS_ratio >= 1/2) &&
            (μS_ratio != 1 || μJ_ratio != 1 || νS_ratio != 1 || νJ_ratio != 1) &&
            !(νJ_ratio == 2 && μS_ratio == 2 && μJ_ratio == 1) && !(νJ_ratio == 2 && μS_ratio == 1 && μJ_ratio == 1/2) &&
            !(μS_ratio == 2 && νJ_ratio == 2 && νS_ratio == 1) && !(μS_ratio == 2 && νJ_ratio == 1 && νS_ratio == 1/2) &&
            !(νJ_ratio == 1/2 && μS_ratio == 1 && μJ_ratio == 2) && !(νJ_ratio == 1/2 && μS_ratio == 1/2 && μJ_ratio == 1) &&
            !(μS_ratio == 1/2 && νJ_ratio == 1 && νS_ratio == 2) && !(μS_ratio == 1/2 && νJ_ratio == 1/2 && νS_ratio == 1)   
        )
            push!(constrained_matrix, (μS_ratio=μS_ratio, μJ_ratio=μJ_ratio, νS_ratio=νS_ratio, νJ_ratio=νJ_ratio))
        end
    end

    return constrained_matrix
end

#Old

#function ratios()
#
#    ratio_matrix = DataFrame(μS_ratio=Float64[], μJ_ratio=Float64[], νS_ratio=Float64[], νJ_ratio=Float64[])
#
#    for i in [1/2, 1, 2]
#        for j in [1/2, 1, 2]
#            for k in [1/2, 1, 2]
#                for l in [1/2, 1, 2]    
#                    push!(ratio_matrix, (μS_ratio=i, μJ_ratio=j, νS_ratio=k, νJ_ratio=l))
#                end
#            end
#        end
#    end
#
#    n_rows = size(ratio_matrix, 1)
#
#    constrained_matrix = DataFrame(μS_ratio=Float64[], μJ_ratio=Float64[], νS_ratio=Float64[], νJ_ratio=Float64[])
#
#    for i in 1:n_rows
#
#        μS_ratio = ratio_matrix[i, :μS_ratio]
#        μJ_ratio = ratio_matrix[i, :μJ_ratio]
#        νS_ratio = ratio_matrix[i, :νS_ratio]
#        νJ_ratio = ratio_matrix[i, :νJ_ratio]
#
#        if (
#            (abs(μS_ratio/μJ_ratio) < 3) && (abs(μS_ratio/νS_ratio) < 3) && (abs(νJ_ratio/νS_ratio) < 3) &&
#            (abs(μS_ratio/μJ_ratio) > 1/3) && (abs(μS_ratio/νS_ratio) > 1/3) && (abs(νJ_ratio/νS_ratio) > 1/3) &&
#            (μS_ratio != 1 || μJ_ratio != 1 || νS_ratio != 1 || νJ_ratio != 1) && 
#            (μJ_ratio != 1/2 || νS_ratio != 1/2 || νJ_ratio != 1) && (μJ_ratio != 1/2 || νJ_ratio != 2 || νS_ratio != 1) 
#        )
#            push!(constrained_matrix, (μS_ratio=μS_ratio, μJ_ratio=μJ_ratio, νS_ratio=νS_ratio, νJ_ratio=νJ_ratio))
#        end
#    end
#
#    return constrained_matrix
#end

function resum_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL)
    ratio_matrix = ratios()
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        μS_ratio = ratio_matrix[i, :μS_ratio]
        μJ_ratio = ratio_matrix[i, :μJ_ratio]
        νS_ratio = ratio_matrix[i, :νS_ratio]
        νJ_ratio = ratio_matrix[i, :νJ_ratio]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL,
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

function fixed_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL)
    ratio_matrix = DataFrame(μS_ratio=Float64[1/2,2], μJ_ratio=Float64[1/2,2], νS_ratio=Float64[1/2,2], νJ_ratio=Float64[1/2,2], μH_ratio=Float64[1/2,2])
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        μS_ratio = ratio_matrix[i, :μS_ratio]
        μJ_ratio = ratio_matrix[i, :μJ_ratio]
        νS_ratio = ratio_matrix[i, :νS_ratio]
        νJ_ratio = ratio_matrix[i, :νJ_ratio]
        μH_ratio = ratio_matrix[i, :μH_ratio]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL,
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

function NP_error(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, XLL::String)

    central_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL)
    ratio_matrix = DataFrame(bmax_ratio=[1/2,2])
    n_rows = size(ratio_matrix, 1)

    delta_list = []

    for i in 1:n_rows
        bmax_ratio = ratio_matrix[i, :bmax_ratio]

        varied_list = sigma_fast(Q=Q, χ_list=χ_list, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, XLL=XLL,
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

