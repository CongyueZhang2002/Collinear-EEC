# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")

# Soft

function γS1_func(nf)

    γS1 = 0

    return γS1
end

function γS2_func(nf) 

    β0 = 4*β0_func(nf) # different convention 1/4π vs 1/π

    γS2 = 2CF * (CA*(-64/9 + 28z3) + β0*(-56/9 + 2z2))

    return γS2
end

function γS3_func(nf)

    β0 = 4*β0_func(nf)
    β1 = 4^2*β0_func(nf)

    γS3 = 2CF * (
          CA^2 * (-37871/162 + 620/27*z2 + 2548/9*z3 + 144z4 - 176/3*z2*z3 - 192z5) 
        + CA*β0 * (4697/54 + 484/27*z2 + 220/9*z3 - 112z4) 
        + β0^2 * (520/81 + 10/3*z2 - 28/3*z3) 
        + β1 * (-1711/54 + 2z2 + 152/9*z3 + 8z4)
    )

    return γS3
end


function γS_func(; αs::Float64, order::Int64, nf::Int64)

    γS1 = γS1_func(nf)
    γS2 = γS2_func(nf)
    γS3 = γS3_func(nf)  

    order1 = αs/(4π) * γS1
    order2 = (αs/(4π))^2 * γS2
    order3 = (αs/(4π))^3 * γS3

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