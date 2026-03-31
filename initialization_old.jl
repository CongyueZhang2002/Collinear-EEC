using Distributed

@everywhere begin

using DelimitedFiles
using CSV
using DataFrames
using Interpolations

#df_Sigma_coeffs = CSV.read("C:\\Users\\congyue zhang\\Desktop\\EEC collinear bspace new\\NP Sigma grid\\Sigma_coeffs_grid.csv", DataFrame)

const a1b_Sigma_vals = unique(df_Sigma_coeffs.a1b)
const a2_Sigma_vals  = unique(df_Sigma_coeffs.a2)

Sigma_coeffs_names = [
    "ΣbLL9",  "ΣbNLL9",  "ΣbNNLL9",
    "ΣbLL8",  "ΣbNLL8",  "ΣbNNLL8",
    "ΣbLL7",  "ΣbNLL7",  "ΣbNNLL7",
    "ΣbLL6",  "ΣbNLL6",  "ΣbNNLL6",
    "ΣbLL5",  "ΣbNLL5",  "ΣbNNLL5",
    "ΣbLL4",  "ΣbNLL4",  "ΣbNNLL4",
    "ΣbLL3",  "ΣbNLL3",  "ΣbNNLL3",
    "ΣbLL2",  "ΣbNLL2",  "ΣbNNLL2",
    "ΣbLL1",  "ΣbNLL1",  "ΣbNNLL1",
    "ΣbLL0",  "ΣbNLL0",  "ΣbNNLL0"
]

Sigma_interpolator = [
    interpolate(
      (a1b_Sigma_vals, a2_Sigma_vals),
      reshape(df_Sigma_coeffs[!, name], length(a2_Sigma_vals), length(a1b_Sigma_vals))',
      Gridded(Linear())
    )
    for name in Sigma_coeffs_names
]

df_Sigma_coeffs = nothing

function Sigma_coeffs_func(; a1b::Float64, a2::Float64) 
    return [f(a1b, a2) for f in Sigma_interpolator]
end

end

#print(Sigma_coeffs_func(a1b=2.3*0.4, a2=1.056))