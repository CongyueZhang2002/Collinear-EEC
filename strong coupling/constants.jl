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

dFAdn = 2.5
dFFdn = 5/36

bdFF = -143.6
bdFA = -998.0
bNCF3 = 80.78
bNCACF2 = -455.247

# Riemann Zeta Function
z2 = 1.64493 
z3 = 1.20206 
z4 = 1.08232
z5 = 1.03693 
z6 = 1.01734
z7 = 1.00835

#
b0 = 1.12292
end

# Renormalon
R0 = 2