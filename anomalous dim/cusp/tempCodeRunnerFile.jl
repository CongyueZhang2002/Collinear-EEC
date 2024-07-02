# cusp anomalous dimension up to 4 loops from https://arxiv.org/pdf/1911.10174.pdf
# expand in order of (αs/4π)^n
using Distributed
include("..\\..\\strong coupling\\constants.jl")
include("..\\..\\strong coupling\\alpha_s.jl")

@everywhere function Γ0_func(nf)
    Γ0 = CF
    return 4*Γ0
end

@everywhere function Γ1_func(nf)
    Γ1 = CF * (1.03864*CA - 0.555556*nf*TF)
    return 4^2*Γ1
end

@everywhere function Γ2_func(nf)
    Γ2 = CF * (1.52982*CA^2 - 1.45614*CA*nf*TF + 0.0562236*CF*nf*TF - 0.0370370*nf^2*TF^2)
    return 4^3*Γ2
end

@everywhere function Γ3_func(nf)
    Γ3 = CF * (
        2.38379*CA^3 - 3.44271*CA^2*nf*TF + 0.303089*CA*CF*nf*TF - 0.242621*CF^2*nf*TF 
        + 0.911990*CA*nf^2*TF^2 - 0.333037*CF*nf^2*TF^2 + 0.0766956*nf^3*TF^3
    ) - dFdAdNA*NA/NF*1.97915 - nf*dFdFdNA*NA/NF*0.483964
    return 4^4*Γ3
end

print(Γ2_func(5))