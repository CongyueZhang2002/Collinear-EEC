using CSV, DataFrames, Interpolations

#path = raw"C:\Users\congyue zhang\Desktop\EEC collinear bspace broken\NP splitting grid\moments_grid.csv"
#df_moments = CSV.read(path, DataFrame)

#Initialize

a1b_vals = unique(df_moments.a1b)
a2_vals  = unique(df_moments.a2)

moments_names = ["γqq0","γqg0","γgq0","γgg0",
              "dγqq0","dγqg0","dγgq0","dγgg0",
              "ddγqq0","ddγqg0","ddγgq0","ddγgg0",
              "γqq1","γqg1","γgq1","γgg1",
              "dγqq1","dγqg1","dγgq1","dγgg1",
              "γqq2","γqg2","γgq2","γgg2"]

interpolator = [
    interpolate(
      (a1b_vals, a2_vals),
      reshape(df_moments[!, name], length(a2_vals), length(a1b_vals))',
      Gridded(Linear())
    )
    for name in moments_names
]

df_moments = nothing

function splitting_moments_func(; b::Float64, a1::Float64, a2::Float64) 
    
    a1b = a1*b

    return [f(a1b, a2) for f in interpolator]
end