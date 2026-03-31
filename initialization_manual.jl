using Distributed

@everywhere begin

using DelimitedFiles

a1b_Sigma_vals = unique(df_Sigma_coeffs[:,1])
a2_Sigma_vals  = unique(df_Sigma_coeffs[:,2])

Na1b, Na2      = 1000, 1000
a1bmin, a1bmax = 0.0, 25.0
a2min, a2max   = 0.5, 2.0
p              = 2

function Sigma_coeffs_func(; a1b::Float64, a2::Float64) 

    x, y = a1b, a2

    i = Int(floor((Na1b-1)*(x/a1bmax)^(1/p)))+1
    j = Int(floor((Na2-1)*((y-a2min)/(a2max-a2min))))+1

    x, y = a1b, a2

    x1, x2 = a1b_Sigma_vals[i], a1b_Sigma_vals[i+1]
    y1, y2 = a2_Sigma_vals[j],  a2_Sigma_vals[j+1]

    row11 = (i-1)*Na2 + j
    row12 = (i-1)*Na2 + (j+1)
    row21 =  i   *Na2 + j
    row22 =  i   *Na2 + (j+1)

    z11 = df_Sigma_coeffs[row11, 3:end]
    z12 = df_Sigma_coeffs[row12, 3:end]
    z21 = df_Sigma_coeffs[row21, 3:end]
    z22 = df_Sigma_coeffs[row22, 3:end]

    tx = (x - x1)/(x2 - x1)
    ty = (y - y1)/(y2 - y1)

    c11 = (1-tx)*(1-ty)
    c12 = (1-tx)*ty
    c21 = tx*(1-ty)
    c22 = tx*ty
    
    result = c11 .* z11 .+ c12 .* z12 .+ c21 .* z21 .+ c22 .* z22

    return result
end

end

#print(Sigma_coeffs_func(a1b=2.3*0.4, a2=1.056))