function hq0_func(nf)

    value = 0.5

    return value
end

function hq1_func(nf)

        value = 10.9167

    return value
end

function hq2_func(nf)

        value = 343.001 - 17.762*nf

    return value
end

function hq_func(; α::Float64, order::Int64, nf::Int64)

    order0 = hq0_func(nf)
    order1 = hq1_func(nf)
    order2 = hq2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    elseif order == 2
        total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end 

#--------------------------------------------------------------------------------------------------

function hg0_func(nf)

        value = 0.0

    return value
end

function hg1_func(nf)

        value = -1.97222

    return value
end

function hg2_func(nf)

        value = -130.028

    return value
end

function hg_func(; α::Float64, order::Int64, nf::Int64)

    order0 = hg0_func(nf)
    order1 = hg1_func(nf)
    order2 = hg2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    elseif order == 2
        total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end 

#--------------------------------------------------------------------------------------------------

function dhq0_func(nf)

        value = 0.0

    return value
end

function dhq1_func(nf)

        value = 2.6255

    return value
end

function dhq_func(; α::Float64, order::Int64, nf::Int64)

    order0 = 0
    order1 = dhq1_func(nf)
    #order2 = dhq2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    #elseif order == 2
    #    total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end 

#--------------------------------------------------------------------------------------------------

function dhg0_func(nf)

        value = 0.0

    return value
end

function dhg1_func(nf)

        value = 1.30394

    return value
end

function dhg_func(; α::Float64, order::Int64, nf::Int64)

    order0 = 0
    order1 = dhg1_func(nf)
    #order2 = dhg2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    #elseif order == 2
    #    total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end

#--------------------------------------------------------------------------------------------------

function ddhq1_func(nf)

    value = -0.762183

    return value
end

function ddhg1_func(nf)

    value = -1.84622

    return value
end

function ddhq_func(; α::Float64, order::Int64, nf::Int64)

    order0 = 0
    order1 = ddhq1_func(nf)
    #order2 = dhg2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    #elseif order == 2
    #    total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end

function ddhg_func(; α::Float64, order::Int64, nf::Int64)

    order0 = 0
    order1 = ddhg1_func(nf)
    #order2 = dhg2_func(nf)

    if order == 0
        total = order0
    elseif order == 1
        total = order0 + (α/(4*π))*order1
    #elseif order == 2
    #    total = order0 + (α/(4*π))*order1 + (α/(4*π))^2*order2
    end

    return total
end