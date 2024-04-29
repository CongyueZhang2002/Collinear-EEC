import Pkg
Pkg.add("HCubature")

using HCubature

# Define the function you want to integrate
f(x) = exp(-x[1]^2)  # x[1] accesses the first element of the array x

# Set the limits of integration as arrays
xmin = [0]  # an array with one element for the lower bound
xmax = [1]  # an array with one element for the upper bound

# Perform the integration
integral, error = hcubature(f, xmin, xmax)

println("Integral: ", integral)
println("Estimated error: ", error)

