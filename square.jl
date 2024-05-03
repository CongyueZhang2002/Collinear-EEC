using Distributed

@everywhere function square(x,b)
    return x^2 + b
end