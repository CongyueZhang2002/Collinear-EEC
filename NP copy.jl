function NP(; b::Float64, Q::Float64, bmax::Float64, parameters::Vector{Float64})

    a = parameters

    if length(a)==0
        NP=1.0
    elseif length(a)==1
        NP = exp(-0.5*a[1]*b^2)
    elseif length(a)==2
        NP = exp(-0.5*a[1]*b^2)*(1-2*a[2]*b)
        #NP = exp(-0.5*a[1]*b^2+a[2]*b) 
    elseif length(a)==3
        NP = exp(-0.5*a[1]*b^2+a[2]*b+a[3]*b^0.5) 
    end  

    bstar = b/(1+(b/bmax)^2)^0.5
    g2 = 0.84
    Q0 = sqrt(2.4)

    return NP#exp(-g2/2*log(b/bstar)*log(Q/Q0))
end

function NP_new(; b::Float64, Q::Float64, μJ::Float64, νdQ::Float64, αs::Float64, bmax::Float64, parameters::Vector{Float64}, order::Int64)

    g = parameters

    if length(g)==0
        NP=1.0
    elseif length(g)==1
        g1=g[1]
    elseif length(g)==2
        g1=g[1]
        g2=g[2]
    elseif length(g)==3
        g1=g[1]
        g2=g[2]
        g3=g[3]
    end  

    a = exp(-g1*b^2)
    a2 = a^2

    log1 = log(a)
    log2 = log1^2
    log3 = log1^3
    log4 = log1^4
    log5 = log1^5
    log6 = log1^6
    log7 = log1^7
    log8 = log1^8
    log9 = log1^9
    log10 = log1^10
    
    A0=(a2*(14175 + 283500*log1 + 1275750*log2 + 2268000*log3 + 1984500*log4 + 952560*log5 + 264600*log6 + 43200*log7 + 
    4050*log8 + 200*log9 + 4*log10))/14175
    A1=(-4*a2*log1*(155925 + 935550*log1 + 1871100*log2 + 1746360*log3 + 873180*log4 + 
    249480*log5 + 41580*log6 + 3960*log7 + 198*log8 + 4*log9))/2835
    A2=(8*a2*log1*(51975 + 363825*log1 + 790020*log2 + 776160*log3 + 401940*log4 + 117810*log5 + 20020*log6 + 1936*log7 + 98*log8 + 
    2*log9))/315
    A3=(-32*a2*log1*(155925 + 1185030*log1 + 2723490*log2 + 2785860*log3 + 1486485*log4 + 445830*log5 + 77154*log6 + 
    7572*log7 + 388*log8 + 8*log9))/945
    A4=(32*a2*log1*(62370 + 498960*log1 + 1193940*log2 + 1260765*log3 + 690030*log4 + 
    211266*log5 + 37188*log6 + 3702*log7 + 192*log8 + 4*log9))/135
    A5= (-128*a2*log1*(51975 + 430650*log1 + 1061775*log2 + 1150050*log3 + 643170*log4 + 200580*log5 + 35870*log6 + 3620*log7 + 
    190*log8 + 4*log9))/225
    A6=(128*a2*log1*(44550 + 378675*log1 + 955350*log2 + 1056330*log3 + 601740*log4 + 190770*log5 + 
    34620*log6 + 3540*log7 + 188*log8 + 4*log9))/135
    A7= (-256*a2*log1*(155925 + 1351350*log1 + 3471930*log2 + 3904740*log3 + 2259810*log4 + 727020*log5 + 133740*log6 + 13848*log7 + 
    744*log8 + 16*log9))/945
    A8=(512*a2*log1*(17325 + 152460*log1 + 397530*log2 + 453495*log3 + 266070*log4 + 86730*log5 + 
    16156*log6 + 1693*log7 + 92*log8 + 2*log9))/315
    A9=(-1024*a2*log1*(31185 + 277830*log1 + 733320*log2 + 846720*log3 + 502740*log4 + 165816*log5 + 31248*log6 + 3312*log7 + 
    182*log8 + 4*log9))/2835
    A10=(2048*a2*log1*(14175 + 127575*log1 + 340200*log2 + 396900*log3 + 238140*log4 + 79380*log5 + 
    15120*log6 + 1620*log7 + 90*log8 + 2*log9))/14175

    Lb = log((b*μJ/b0)^2)
    LQ = log(νdQ)

    J0 = a
    J1 = αs/(4π)*(
        8*(5*A1 + (2429702990*A10 + 363*(5953500*A2 + 6228096*A3 + 6379800*A4 + 6477144*A5 + 6545469*A6 + 6596344*A7 + 6635838*A8) + 2420288522*A9)/384199200 + 
        ((83160*A1 + 169312*A10 + 33*(3150*A2 + 3584*A3 + 3920*A4 + 4196*A5 + 4431*A6 + 4636*A7 + 4818*A8) + 164398*A9)*Lb)/27720 - 
        (2*(A1 + A10 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9)*π^2)/3 + (A0*(12 + 9*Lb - 4*π^2))/6 - (A0 + A1 + A10 + A2 + A3 + A4 + A5 + A6 + A7 + A8 + A9)*Lb*LQ)
    )/3

    LΛ = log(0.25^2*b^2/b0^2)

    C = -A0 - (8*A1)/3 - (7591*A10)/1260 - (7*A2)/2 - (61*A3)/15 - (9*A4)/2 - (1019*A5)/210 - (103*A6)/20 - (3407*A7)/630 - (789*A8)/140 - (80939*A9)/13860 - 
        2/3 - LΛ 

    if order == 1 
        total = J0 + J1 + a*g2*b^2#*C
    elseif order == 2  
        total = J0 + J1 + J2 + a*g2*b^2#*C
    elseif order == 3 
        total = J0 + J1 + J2 + J3 + a*g2*b^2#*C
    end

    return total
end