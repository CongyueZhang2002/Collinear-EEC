using Distributed
using HCubature
using SpecialFunctions

@everywhere begin
    
    using HCubature
    using SpecialFunctions

    include("strong coupling\\constants.jl")
    include("strong coupling\\alpha_s.jl")
    include("anomalous dim\\cusp\\cusp.jl")
    include("anomalous dim\\non-cusp\\gammaH.jl")
    include("anomalous dim\\non-cusp\\gammaB.jl")
    include("anomalous dim\\non-cusp\\gammaS.jl")
    include("anomalous dim\\non-cusp\\gammav.jl")
    include("hard\\hard.jl")
    include("jet\\jet.jl")
    include("soft\\soft.jl")

    function Integral_H(; Q::Float64, μB::Float64, μH::Float64, 
                        αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

        f(x)=1/x[1]*γH_final(Q=Q, μ=x[1], αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
        integral, error = hcubature(f, [μH], [μB])

        return integral
    end

    #function Integral_J(; Q::Float64, νJ::Float64, μ::Float64, μJ::Float64, 
    #                    αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)
    #
    #    f(x)=1/x[1]*γJ_final(μ=x[1], νdQ=νJ/Q, αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
    #    integral, error = hcubature(f, [μJ], [μ])
    #
    #    return integral
    #end

    function Integral_S(; νS::Float64, μB::Float64, μS::Float64, 
                        αs_ini::Float64, μ_ini::Float64, order::Int64, nf::Int64)

        f(x)=1/x[1]*γS_final(μ=x[1], ν=νS, αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)
        integral, error = hcubature(f, [μS], [μB])
        return integral
    end

    function Integrand(;b::Float64, Q::Float64, z::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
        μJ_ratio::Float64, νJ_ratio::Float64, μH_ratio::Float64,
        μS_ratio::Float64, νS_ratio::Float64, bmax_ratio::Float64)

        if b < 1.0e-9
            return 0.0
        end

        J0 = besselj0(b*Q*(1-z)^0.5)

        bmax = bmax_ratio * b0
        bstar = b/(1+(b/bmax)^2)^0.5

        μJ = μJ_ratio * b0/bstar
        νJ = νJ_ratio * Q
        μS = μS_ratio * b0/bstar
        νS = νS_ratio * b0/b
        μH = μH_ratio * Q    

        αs_J = alpha_s_func(μf=μJ, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
        αs_S = alpha_s_func(μf=μS, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)

        J = J_func(b=b, μJ=μJ, νdQ=νJ/Q, αs=αs_J, order=order, nf=nf)
        S = S_func(b=b, μS=μS, νS=νS, αs=αs_S, order=order, nf=nf)

        γH_Integral = Integral_H(Q=Q, μB=μJ, μH=μH, 
                                αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)

        γS_Integral = Integral_S(νS=νS, μB=μJ, μS=μS, 
                                αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf)

        γν_Integral = γν_final(b=b, μB=μJ, 
                                αs_ini=αs_ini, μ_ini=μ_ini, order=order, nf=nf, bmax=bmax)

        total = Q^2/4*b*J0*J*J*S*exp(γH_Integral+γS_Integral)*(νJ/νS)^γν_Integral

        return total
    end

    function Hard_Part(; Q::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, μH_ratio::Float64, order::Int64)
        
        μH = μH_ratio * Q
        αs_H = alpha_s_func(μf=μH, μi=μ_ini, αs=αs_ini, order=order+1, nf=nf)
        H = H_func(μH=μH, Q=Q, αs=αs_H, order=order)

        return H
    end

    function sigma_z(; Q::Float64, z::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
                    μH_ratio::Float64, μJ_ratio::Float64, νJ_ratio::Float64, 
                    μS_ratio::Float64, νS_ratio::Float64, bmax_ratio::Float64)

        f(x) = Integrand(b=x[1], Q=Q, z=z, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                        μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, μS_ratio=μS_ratio, 
                        νS_ratio=νS_ratio, μH_ratio=μH_ratio, bmax_ratio=bmax_ratio)

        integral, error = hcubature(f, [0.0], [10.0])
        
        hard = Hard_Part(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, μH_ratio=μH_ratio, order=order)

        total = hard*integral

        return total
    end

    function sigma_χ(; Q::Float64, χ::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
        μH_ratio::Float64=1.0, μJ_ratio::Float64=1.0, νJ_ratio::Float64=1.0, 
        μS_ratio::Float64=1.0, νS_ratio::Float64=1.0, bmax_ratio::Float64=1.0)

        z = 0.5*(1-cos(χ))

        part = sigma_z(Q=Q, z=z, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                        μH_ratio=μH_ratio, μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, 
                        μS_ratio=μS_ratio, νS_ratio=νS_ratio, bmax_ratio=bmax_ratio)

        total = 0.5*sin(χ)*part

        return total
    end

end

function sigma_fast(; Q::Float64, χ_list::Vector{Float64}, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
    μH_ratio::Float64=1.0, μJ_ratio::Float64=1.0, νJ_ratio::Float64=1.0, 
    μS_ratio::Float64=1.0, νS_ratio::Float64=1.0, bmax_ratio::Float64=1.0)
    try
        #if length(workers()) < cores
        #    addprocs(cores - length(workers()))  
        #end

        Q1 = Q
        μ_ini1 = μ_ini
        αs_ini1 = αs_ini
        nf1 = nf
        order1 = order

        function compute_sigma_χ(χ)
            return sigma_χ(Q=Q1, χ=χ, μ_ini=μ_ini1, αs_ini=αs_ini1, nf=nf1, order=order1,
                            μH_ratio=μH_ratio, μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, 
                            μS_ratio=μS_ratio, νS_ratio=νS_ratio, bmax_ratio=bmax_ratio)
        end

        results = pmap(compute_sigma_χ, χ_list)

        return results
    catch e
        println("An error occurred: ", e)
        rethrow(e)
    finally
        #rmprocs(workers()) 
    end
end

# For checking dσ/dz=∱db (integrated)
function integrated(; Q::Float64, b::Float64, μ_ini::Float64, αs_ini::Float64, nf::Int64, order::Int64,
    μH_ratio::Float64 = 1.0, μJ_ratio::Float64 = 1.0, νJ_ratio::Float64 = 1.0, 
    μS_ratio::Float64 = 1.0, νS_ratio::Float64 = 1.0, bmax_ratio::Float64 = 1.0)

    f = Integrand(b=b, Q=Q, z=1.0, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, order=order,
                        μJ_ratio=μJ_ratio, νJ_ratio=νJ_ratio, μH_ratio=μH_ratio, 
                        μS_ratio=μS_ratio, νS_ratio=νS_ratio, bmax_ratio=bmax_ratio)

    hard = Hard_Part(Q=Q, μ_ini=μ_ini, αs_ini=αs_ini, nf=nf, μH_ratio=μH_ratio, order=order)

    total = hard*f

    return total
end