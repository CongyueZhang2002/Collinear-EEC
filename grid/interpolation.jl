using Interpolations
using DelimitedFiles
using HCubature

include("..\\strong coupling\\constants.jl")
include("..\\strong coupling\\alpha_s.jl")
include("..\\anomalous dim\\cusp\\cusp.jl")

@everywhere begin
    μf_array = range(log(0.50),log(2500.0),5000)
    αs_ini_array = 0.1129:0.00005:0.1231

    #Γoverμ_1(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.2, order=1, nf=5)
    #Γoverμ_2(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.2, order=2, nf=5)
    #Γoverμ_3(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.2, order=3, nf=5)
    #Γoverμ_4(μf,αs_ini) = Γoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.2, order=4, nf=5)

    #Grid1 = readdlm("Grid/Γoverμ_1.csv", ',')
    #Grid2 = readdlm("Grid/Γoverμ_2.csv", ',')
    #Grid3 = readdlm("Grid/Γoverμ_3.csv", ',')
    #Grid4 = readdlm("Grid/Γoverμ_4.csv", ',')

    #print(Γoverμ_1(1.125,0.118025),interp_cubic_1(1.125, 0.118025)-Γoverμ_1(1.125,0.118025))
    #print(Γoverμ_2(1.125,0.118025),interp_cubic_2(1.125, 0.118025)-Γoverμ_2(1.125,0.118025))
    #print(Γoverμ_3(1.125,0.118025),interp_cubic_3(1.125, 0.118025)-Γoverμ_3(1.125,0.118025))
    #print(Γoverμ_4(1.125,0.118025),interp_cubic_4(1.125, 0.118025)-Γoverμ_4(1.125,0.118025))

    interp_cubic_1 = cubic_spline_interpolation((μf_array, αs_ini_array), Grid1)
    interp_cubic_2 = cubic_spline_interpolation((μf_array, αs_ini_array), Grid2)
    interp_cubic_3 = cubic_spline_interpolation((μf_array, αs_ini_array), Grid3)
    interp_cubic_4 = cubic_spline_interpolation((μf_array, αs_ini_array), Grid4)

    function Γ_integral(; μi::Float64, μf::Float64, αs_ini::Float64, order::Int64)

        μ_max = μf
        μ_min = μi

        if order == 1 
            inter = interp_cubic_1
        elseif order == 2  
            inter = interp_cubic_2
        elseif order == 3 
            inter = interp_cubic_3
        elseif order == 4 
            inter = interp_cubic_4
        end

        return inter(log(μf),αs_ini)-inter(log(μi),αs_ini)

    end

    #print(Γoverμ_integral(μf=91.2, αs_ini=0.118, μ_ini=91.1876, order=4, nf=5)-Γoverμ_integral(μf=1.0, αs_ini=0.118, μ_ini=91.1876, order=4, nf=5))
    #print(Γ_integral(μi=1.0, μf=91.2, αs_ini=0.118, order=4))

    #---------------------------------------------------------------------------------------

    #lnΓoverμ_1(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=1, nf=5)
    #lnΓoverμ_2(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=2, nf=5)
    #lnΓoverμ_3(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=3, nf=5)
    #lnΓoverμ_4(μf,αs_ini) = lnΓoverμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, order=4, nf=5)

    #lnGrid1 = readdlm("Grid/lnΓoverμ_1.csv", ',')
    #lnGrid2 = readdlm("Grid/lnΓoverμ_2.csv", ',')
    #lnGrid3 = readdlm("Grid/lnΓoverμ_3.csv", ',')
    #lnGrid4 = readdlm("Grid/lnΓoverμ_4.csv", ',')

    lninterp_cubic_1 = cubic_spline_interpolation((μf_array, αs_ini_array), lnGrid1)
    lninterp_cubic_2 = cubic_spline_interpolation((μf_array, αs_ini_array), lnGrid2)
    lninterp_cubic_3 = cubic_spline_interpolation((μf_array, αs_ini_array), lnGrid3)
    lninterp_cubic_4 = cubic_spline_interpolation((μf_array, αs_ini_array), lnGrid4)

    function lnΓ_integral(; μi::Float64, μf::Float64, αs_ini::Float64, order::Int64)

        μ_max = μf
        μ_min = μi

        if order == 1 
            inter = lninterp_cubic_1
        elseif order == 2  
            inter = lninterp_cubic_2
        elseif order == 3 
            inter = lninterp_cubic_3
        elseif order == 4 
            inter = lninterp_cubic_4
        end

        return inter(log(μf),αs_ini)-inter(log(μi),αs_ini)

    end

    #---------------------------------------------------------------------------------------

    #α1overμ(μf,αs_ini) = α1overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
    #α2overμ(μf,αs_ini) = α2overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
    #α3overμ(μf,αs_ini) = α3overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)
    #α4overμ(μf,αs_ini) = α4overμ_integral(μf=μf, αs_ini=αs_ini, μ_ini=91.1876, nf=5)

    #α1Grid = readdlm("Grid/α1overμ.csv", ',')
    #α2Grid = readdlm("Grid/α2overμ.csv", ',')
    #α3Grid = readdlm("Grid/α3overμ.csv", ',')
    #α4Grid = readdlm("Grid/α4overμ.csv", ',')

    α1interp_cubic = cubic_spline_interpolation((μf_array, αs_ini_array), α1Grid)
    α2interp_cubic = cubic_spline_interpolation((μf_array, αs_ini_array), α2Grid)
    α3interp_cubic = cubic_spline_interpolation((μf_array, αs_ini_array), α3Grid)
    α4interp_cubic = cubic_spline_interpolation((μf_array, αs_ini_array), α4Grid)

    function αn_integral(; μi::Float64, μf::Float64, αs_ini::Float64, power::Int64)

        μ_max = μf
        μ_min = μi

        if power == 1 
            inter = α1interp_cubic
        elseif power == 2  
            inter = α2interp_cubic
        elseif power == 3 
            inter = α3interp_cubic
        elseif power == 4 
            inter = α4interp_cubic
        end

        return inter(log(μf),αs_ini)-inter(log(μi),αs_ini)

    end
end
#print(α4overμ_integral(μf=1.0026, αs_ini=0.118, μ_ini=91.1876,nf=5)-αn_integral(μi=1.0, μf=1.0026, αs_ini=0.118, power=4))
#print(αn_integral(μi=1.0, μf=1.0026, αs_ini=0.118, power=4))