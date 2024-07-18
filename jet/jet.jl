# jet up to 3 loops from https://arxiv.org/abs/2012.07859

using Distributed
include("..\\strong coupling\\constants.jl")

# jet function

@everywhere function J1_func(αs,Lb,LQ,nf)

    J1 = (
        1 + (CF*αs)/π + (3*CF*Lb*αs)/(4*π) + 
        (CF*Lb*LQ*αs)/π - (CF*π*αs)/3
    )

    return J1
end

@everywhere function J2_func(αs,Lb,LQ,nf)

    J2 = αs^2*(
        -0.52108 - 0.6645800000000001*Lb + 0.14776*Lb^2 - 0.26965*LQ - 0.17915*Lb*LQ + 
        0.26456*Lb^2*LQ + 0.09006*Lb^2*LQ^2      
    )

    return J2
end

@everywhere function J3_func(αs,Lb,LQ,nf)

    J3 = αs^3*(
        -0.42938707833771267 - 1.0417604864269598*Lb - 0.434427523774156*Lb^2 + 
        0.07577617465207752*Lb^3 - 0.08267770782339488*LQ - 0.7422117758365783*Lb*LQ - 
        0.13930697059087566*Lb^2*LQ + 0.1565792088377247*Lb^3*LQ - 0.11444278880208449*Lb*LQ^2 + 
        0.011493849910393553*Lb^2*LQ^2 + 0.08361508927125795*Lb^3*LQ^2 + 
        0.01274134693657264*Lb^3*LQ^3
    )

    return J3
end

@everywhere function J_func(; b::Float64, μJ::Float64, νdQ::Float64, αs::Float64, order::Int64, nf::Int64)

    Lb = log((b*μJ/b0)^2)
    LQ = log(νdQ)

    J1 = J1_func(αs,Lb,LQ,nf)
    J2 = J2_func(αs,Lb,LQ,nf)
    J3 = J3_func(αs,Lb,LQ,nf)

    if order == 1 
        total = J1
    elseif order == 2  
        total = J1 + J2
    elseif order == 3 
        total = J1 + J2 + J3
    end

    return total

end
