function NP(; b::Float64, parameters::Vector{Float64})

    a = parameters

    # Choose your parameterization 

    NP = 1.0
    NP = exp(-0.5*a[1]*b^2)
    #NP = exp(-0.5*a1*b^2)*(1-2*a2*b)
    #NP = exp(-0.5*a1*b^2+a2*b+a3*b^0.5)    
    print(a[0])
    return NP
end

NP(b=1.0, parameters=[1.0])