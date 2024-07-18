function NP(; b::Float64, parameters::Vector{Float64})

    a = parameters

    if length(a)==0
        NP=1.0
    elseif length(a)==1
        NP = exp(-0.5*a[1]*b^2)
    elseif length(a)==2
        NP = exp(-0.5*a[1]*b^2)*(1-2*a[2]*b)
    elseif length(a)==3
        NP = exp(-0.5*a[1]*b^2+a[2]*b+a[3]*b^0.5)
    end  
    return NP
end
