using Distributed

#@everywhere begin

using Interpolations

# load splitting grid

a1b_splitting_vals = unique(df_splitting_moments.a1b)
a2_splitting_vals  = unique(df_splitting_moments.a2)

splitting_moments_names = ["γqq0","γqg0","γgq0","γgg0",
                        "dγqq0","dγqg0","dγgq0","dγgg0",
                        "ddγqq0","ddγqg0","ddγgq0","ddγgg0",
                        "γqq1","γqg1","γgq1","γgg1",
                        "dγqq1","dγqg1","dγgq1","dγgg1",
                        "γqq2","γqg2","γgq2","γgg2"]

splitting_interpolator = [
    interpolate(
      (a1b_splitting_vals, a2_splitting_vals),
      reshape(df_splitting_moments[!, name], length(a2_splitting_vals), length(a1b_splitting_vals))',
      Gridded(Linear())
    )
    for name in splitting_moments_names
]

df_splitting_moments = nothing

# load hard grid

a1b_hard_vals = unique(df_hard_moments.a1b)
a2_hard_vals  = unique(df_hard_moments.a2)

hard_moments_names = ["hq0","hq1","hq2",
                    "hg0","hg1","hg2"]

hard_interpolator = [
    interpolate(
      (a1b_hard_vals, a2_hard_vals),
      reshape(df_hard_moments[!, name], length(a2_hard_vals), length(a1b_hard_vals))',
      Gridded(Linear())
    )
    for name in hard_moments_names
]

df_hard_moments = nothing

function splitting_moments_func(; a1b::Float64, a2::Float64) 
    return [f(a1b, a2) for f in splitting_interpolator]
end

function hard_moments_func(; a1b::Float64, a2::Float64) 
    return [f(a1b, a2) for f in hard_interpolator]
end

function init_splitting_globals(; a1b::Float64, a2::Float64)
    vals = splitting_moments_func(a1b=a1b, a2=a2)

    global γqq0NP, γqg0NP, γgq0NP, γgg0NP,
           dγqq0NP, dγqg0NP, dγgq0NP, dγgg0NP,
           ddγqq0NP, ddγqg0NP, ddγgq0NP, ddγgg0NP,
           γqq1NP, γqg1NP, γgq1NP, γgg1NP,
           dγqq1NP, dγqg1NP, dγgq1NP, dγgg1NP,
           γqq2NP, γqg2NP, γgq2NP, γgg2NP

    (γqq0NP, γqg0NP, γgq0NP, γgg0NP,
     dγqq0NP, dγqg0NP, dγgq0NP, dγgg0NP,
     ddγqq0NP, ddγqg0NP, ddγgq0NP, ddγgg0NP,
     γqq1NP, γqg1NP, γgq1NP, γgg1NP,
     dγqq1NP, dγqg1NP, dγgq1NP, dγgg1NP,
     γqq2NP, γqg2NP, γgq2NP, γgg2NP) = vals

    return nothing
end

function init_hard_globals(; a1b::Float64, a2::Float64)
    vals = hard_moments_func(a1b=a1b, a2=a2)

    global hq0NP,hq1NP,hq2NP,hg0NP,hg1NP,hg2NP

    (hq0NP,hq1NP,hq2NP,hg0NP,hg1NP,hg2NP) = vals

    return nothing
end

#end