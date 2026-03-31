using CSV, DataFrames, ProgressMeter

# ── 1) Grid parameters
Na1b, Na2       = 1000, 1000
a1bmin, a1bmax = 0.0, 25.0
a2min, a2max   = 0.5, 2.0
p               = 2

a1b_vals = a1bmin .+ (a1bmax - a1bmin) .* collect(range(0, 1; length=Na1b)).^p
# Linear spacing for a2
a2_vals  = collect(range(a2min, a2max; length=Na2))

#y = 24.949975000025052*10^(-5)
#print(Int(floor((Na1b-1)*(y/a1bmax)^(1/p)))+1)
#print(Int(floor((Na1b-1)*(y/a1bmax)^(1/p)))+1)

#y=1.9985
#print(Int(floor((Na2-1)*((y-a2min)/(a2max-a2min))))+1)

function billinear_Sigma(; a1b::Float64, a2::Float64)

    Na1b, Na2      = 1000, 1000
    a1bmin, a1bmax = 0.0, 25.0
    a2min, a2max   = 0.5, 2.0
    p              = 2

    x, y = a1b, a2

    i = Int(floor((Na1b-1)*(y/a1bmax)^(1/p)))+1
    j = Int(floor((Na2-1)*((y-a2min)/(a2max-a2min))))+1

    x, y = a1b, a2

    x1, x2 = a1b_vals[i], a1b_vals[i+1]
    y1, y2 = a2_vals[j],  a2_vals[j+1]

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

#using DelimitedFiles
#df_Sigma_coeffs, header = readdlm("C:\\Users\\congyue zhang\\Desktop\\EEC collinear bspace new\\NP Sigma grid\\Sigma_coeffs_grid.csv", ',', Float64, header=true)
#display(df_Sigma_coeffs[1, 3:end])