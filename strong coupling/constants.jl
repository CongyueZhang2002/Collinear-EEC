using Distributed

@everywhere begin
# Color Factors
CA = 3
CF = 4/3
NC = 3

NA = 8
NF = 3
TF = 1/2

dFdAdNA = 15/16
dFdFdNA = 5/96

# Riemann Zeta Function
z2 = 1.64493 
z3 = 1.20206 
z4 = 1.08232
z5 = 1.03693 
z6 = 1.01734

#
b0 = 1.12292
end

# Renormalon
R0 = 2