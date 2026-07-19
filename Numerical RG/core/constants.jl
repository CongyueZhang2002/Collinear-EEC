# Physical constants adapted from the WSL TMD fitter Core/constants.jl.
#
# The fitter obtains these masses from its configured PDF set. Numerical RG is
# standalone, so these temporary values can be replaced here without changing
# the solver API.
const mc = 1.27
const mb = 4.18
const mt = 172.76

# TMD
const b0 = 1.12292

# SU(3)
const NC = 3
const CF = 4.0 / 3.0
const CA = 3
const TF = 1.0 / 2.0

const NA = NC^2 - 1
const NF = NC # Fundamental-representation dimension, not active-flavor nf.

# SU(3) Casimir tensors
const dAA_NA = 135.0 / 8.0
const dFA_NA = 15.0 / 16.0
const dFF_NA = 5.0 / 96.0

const dFA_NF = dFA_NA * NA / NF
const dFF_NF = dFF_NA * NA / NF

# SU(3) four-loop colour coefficients from arXiv:2205.02249, Table 1
const b_d4_FA = 998.0
const b_d4_FF = 143.6
const b_CA_CF2_nf = 455.247
const b_CF3_nf = -80.780

# Strong-coupling reference values
const MZ = 91.1876
const DEFAULT_ALPHA_S_MZ = 0.118
const ΛQCD3 = 0.326
const ΛQCD4 = 0.326
const ΛQCD5 = 0.226

# Quark charges and electromagnetic coupling
const eu = 2.0 / 3.0
const ed = 1.0 / 3.0
const αem = 1.0 / 137.0

# Riemann zeta constants
const z2 = 1.644934067
const z3 = 1.202056903
const z4 = 1.082323234
const z5 = 1.036927755
const z6 = 1.017343062
const z7 = 1.008349277
