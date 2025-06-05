using SpecialFunctions

a = Dict{Int,Float64}()
a[1]=1.0
a[2]=0.5
b=1.0
print((a[2]*b)^(-a[1])*besselj(a[1], a[2]*b)/((2^(-a[1]))/gamma(1 + a[1])))
print("  ")
print((a[2]*b)^a[1]*besselk(a[1], a[2]*b)/(2^(a[1]-1)*gamma(a[1])))