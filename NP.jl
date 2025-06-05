@everywhere begin

function NP(; z::Float64, b::Float64, ζi::Float64, ζf::Float64, params::Vector{Float64})

    a = params
    bstar = b/(1+(b/b0)^2)^0.5

    if length(a)==0
        NP = 1.0
    elseif length(a)==1
        NP = exp(-a[1]*b)
    elseif length(a)==2
        NP = exp(-a[1]*b^a[2])
        #NP = exp(-a[1]*b-a[2]*b^2)
        #NP = exp(-a[1]*b)*(1+a[2]*b)
        #NP = exp(-(a[1]*b+a[2]))
    elseif length(a)==3
        #NP = exp(-a[1]*b^a[2])*exp(-a[3]*b*bstar*log(ζf/ζi))
        #NP = exp(-(a[1]*b^2+a[2]*b+a[3]*sqrt(b)))
        #NP = exp(-(a[1]*b^2+a[2]*b+a[3]))
        #NP = exp(-(a[1]*b^a[2]+a[3]))
        NP = a[1]*exp(-a[2]*b^a[3])
    elseif length(a)==4

        Ω1 = a[1]
        z0 = a[2]
        p = a[3]
        Ω0 = a[4]

        NP = exp(-(
                    Ω1/(1+(z/z0)^p) + Ω0
                )*b)
        NP = exp(-(exp(-a[3]*z)*a[1]+(1-exp(-a[3]*z))*a[2])*b^a[4])
    end

    return NP
end 

function NP_qg(; b::Float64, ζi::Float64, ζf::Float64, params::Vector{Float64})

    a = params
    bstar = b/(1+(b/b0)^2)^0.5

    if length(a)==0
        NP = 1.0
    elseif length(a)==1
        NP = exp(-a[1]*b)
    elseif length(a)==2
        NP_quark = exp(-a[1]*b)
        NP_gluon = exp(-a[2]*b)
    elseif length(a)==4
        #NP_quark = exp(-a[1]*b^a[3])
        #NP_gluon = exp(-a[2]*b^a[4])
        NP_quark = exp(-a[1]*b)*exp(-a[3]*b*bstar*log(ζf/ζi))
        NP_gluon = exp(-a[2]*b)*exp(-a[4]*b*bstar*log(ζf/ζi))        
    end

    return NP_quark, NP_gluon
end 

end