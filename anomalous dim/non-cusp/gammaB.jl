# non-cusp anomalous dimension up to 3 loops from https://arxiv.org/abs/1909.00811

include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")
include("..\\cusp\\cusp.jl")
include("gammaS.jl")

# Beam

function γB1_func(nf)

    γB1 = 6CF

    return γB1
end

function γB2_func(nf) 

    β0 = 4*β0_func(nf) # different convention 1/4π vs 1/π

    γB2 = 2CF * (
          CA * (73/9 - 40z3) 
        + CF * (3/2 - 12z2 + 24z3) 
        + β0 * (121/18 + 2z2)
    )

    return γB2
end

function γB3_func(nf)

    β0 = 4*β0_func(nf)
    β1 = 4^2*β1_func(nf)

    γB3 = 2CF * (
          CA^2 * (52019/162 - 1682/27*z2 - 2056/9*z3 - 820/3*z4 + 176/3*z2*z3 + 232z5) 
        + CA*CF * (151/4 - 410/3*z2 + 844/3*z3 - 494/3*z4 + 16*z2*z3 + 120z5) 
        + CF^2 * (29/2 + 18z2 + 68z3 + 144z4 - 32z2*z3 - 240z5) 
        + CA*β0 * (-7739/54 + 650/27*z2 - 1276/9*z3 + 617/3*z4) 
        + β0^2 * (-3457/324 + 10/3*z2 + 16/3*z3) 
        + β1 * (1166/27 - 16/3*z2 + 52/9*z3 - 82/3*z4)
    )

    return γB3
end


function γB_func(; αs::Float64, order::Int64, nf::Int64)

    γB1 = γB1_func(nf)
    γB2 = γB2_func(nf)
    γB3 = γB3_func(nf)  

    order1 = αs/(4π) * γB1
    order2 = (αs/(4π))^2 * γB2
    order3 = (αs/(4π))^3 * γB3

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

# final γJ in integration
function γJ_final(; μ::Float64, νdQ::Float64, αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

    αs = alpha_s_func(μf=μ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
    Γ = Γ_func(αs=αs, order=order+1, nf=nf)
    γJ = γB_func(αs=αs, order=order, nf=nf) + γS_func(αs=αs, order=order, nf=nf)

    γJ_f = 2*Γ*log(νdQ) + γJ

    return γJ_f
end