# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")

# rapidity ν

function γν1_func(nf)

    γν1 = 0

    return γν1
end

function γν2_func(nf) 

    β0 = 4*β0_func(nf) # different convention 1/4π vs 1/π

    γν2 = 2CF * (CA*(-64/9 + 28z3) - β0*56/9)

    return γν2
end

function γν3_func(nf)

    β0 = 4*β0_func(nf)
    β1 = 4^2*β0_func(nf)

    γν3 = 2CF * (
        + CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
        + CA*β0 * (3865/54 + 412/27*z2 + 220/9*z3 - 50z4) 
        + β0^2 * (-464/81 - 8z3) 
        + β1 * (-1711/54 + 152/9*z3 + 8z4)
    )

    return γν3
end


function γν_func(; αs::Float64, order::Int64, nf::Int64)

    γν1 = γν1_func(nf)
    γν2 = γν2_func(nf)
    γν3 = γν3_func(nf)  

    order1 = αs/(4π) * γν1
    order2 = (αs/(4π))^2 * γν2
    order3 = (αs/(4π))^3 * γν3

    if order == 1 
        total = order1
    end
    if order == 2  
        total = order1 + order2
    end
    if order == 3 
        total = order1 + order2 + order3
    end

    return total

end