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

print(ratios())