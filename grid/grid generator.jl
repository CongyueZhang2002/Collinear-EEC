using Interpolations
using DelimitedFiles
using HCubature

include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")
include("..\\anomalous dim\\cusp\\cusp.jl")

#Γoverμ_1(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=1, nf=5)
#Γoverμ_2(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=2, nf=5)
#Γoverμ_3(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=3, nf=5)
#Γoverμ_4(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=4, nf=5)

#μf_array = range(log(0.50),log(2500.0),5000)
#αs_ini_array = 0.1129:0.00005:0.1231

#A1 = [Γoverμ_1(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#A2 = [Γoverμ_2(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#A3 = [Γoverμ_3(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#A4 = [Γoverμ_4(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#writedlm("Grid/Γoverμ_1.csv", A1, ',')
#writedlm("Grid/Γoverμ_2.csv", A2, ',')
#writedlm("Grid/Γoverμ_3.csv", A3, ',')
#writedlm("Grid/Γoverμ_4.csv", A4, ',')

#print("1")

#-------------------------------------------------------------------------------------------

#lnΓoverμ_1(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=1, nf=5)
#lnΓoverμ_2(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=2, nf=5)
#lnΓoverμ_3(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=3, nf=5)
#lnΓoverμ_4(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=4, nf=5)

#μf_array = range(log(0.50),log(2500.0),5000)
#αs_ini_array = 0.1129:0.00005:0.1231

#B1 = [lnΓoverμ_1(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#B2 = [lnΓoverμ_2(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#B3 = [lnΓoverμ_3(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#B4 = [lnΓoverμ_4(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
#writedlm("Grid/lnΓoverμ_1.csv", B1, ',')
#writedlm("Grid/lnΓoverμ_2.csv", B2, ',')
#writedlm("Grid/lnΓoverμ_3.csv", B3, ',')
#writedlm("Grid/lnΓoverμ_4.csv", B4, ',')

#print("2")

@everywhere function α1overμ_integral(; μf::Float64, αs_ini::Float64, μ_ini::Float64, nf::Int64)

    μ_max = [μf]
    μ_min = [1.0] 

    f(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)/x[1]

    integral, error = hcubature(f, μ_min, μ_max)

    return integral

end

@everywhere function α2overμ_integral(; μf::Float64, αs_ini::Float64, μ_ini::Float64, nf::Int64)

    μ_max = [μf]
    μ_min = [1.0] 

    f(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)^2/x[1]

    integral, error = hcubature(f, μ_min, μ_max)

    return integral

end

@everywhere function α3overμ_integral(; μf::Float64, αs_ini::Float64, μ_ini::Float64, nf::Int64)

    μ_max = [μf]
    μ_min = [1.0] 

    f(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)^3/x[1]

    integral, error = hcubature(f, μ_min, μ_max)

    return integral

end

@everywhere function α4overμ_integral(; μf::Float64, αs_ini::Float64, μ_ini::Float64, nf::Int64)

    μ_max = [μf]
    μ_min = [1.0] 

    f(x) = alpha_s_func(μf=x[1], μi=μ_ini, αs=αs_ini, order=4, nf=nf)^4/x[1]

    integral, error = hcubature(f, μ_min, μ_max)

    return integral

end

α1overμ(μf,αs_ini) = α1overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
α2overμ(μf,αs_ini) = α2overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
α3overμ(μf,αs_ini) = α3overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
α4overμ(μf,αs_ini) = α4overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)

μf_array = range(log(0.50),log(2500.0),5000)
αs_ini_array = 0.1129:0.00005:0.1231

C1 = [α1overμ(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
C2 = [α2overμ(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
C3 = [α3overμ(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
C4 = [α4overμ(exp(μf),αs_ini) for μf in μf_array, αs_ini in αs_ini_array]
writedlm("Grid/α1overμ.csv", C1, ',')
writedlm("Grid/α2overμ.csv", C2, ',')
writedlm("Grid/α3overμ.csv", C3, ',')
writedlm("Grid/α4overμ.csv", C4, ',')

print("3")