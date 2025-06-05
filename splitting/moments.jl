#qq_x^2--------------------------------------------------------------------------------------------

function γqq0_func(nf)

    if nf == 5
        value = 5.55555
    end

    return value
end

function γqq1_func(nf)

    if nf == 5
        value = 32.4385 
    end

    return value
end

function γqq2_func(nf)

    if nf == 5
        value = 50.2731 #50.2955 maximum error < 0.05%
    end

    return value
end

#qg_x^2--------------------------------------------------------------------------------------------

function γqg0_func(nf)

    if nf == 5
        value = -2.33333
    end

    return value
end

function γqg1_func(nf)

    if nf == 5
        value = 19.7871 
    end

    return value
end

function γqg2_func(nf)

    if nf == 5
        value = -130.585 #-130.6
    end

    return value
end

#gq_x^2--------------------------------------------------------------------------------------------

function γgq0_func(nf)

    if nf == 5
        value = -1.55555
    end

    return value
end

function γgq1_func(nf)

    if nf == 5
        value = -44.8855
    end

    return value
end

function γgq2_func(nf)

    if nf == 5
        value = -373.571 # -373.553
    end

    return value
end

#gg_x^2--------------------------------------------------------------------------------------------

function γgg0_func(nf)

    if nf == 5
        value = 11.7333
    end

    return value
end

function γgg1_func(nf)

    if nf == 5
        value = -79.9426
    end

    return value
end

function γgg2_func(nf)

    if nf == 5
        value = -675.768 # -675.797
    end

    return value
end

function γ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = γqq0_func(nf)
        order1 = γqq1_func(nf)
        order2 = γqq2_func(nf)
    elseif type == "qg"
        order0 = γqg0_func(nf)
        order1 = γqg1_func(nf)
        order2 = γqg2_func(nf)
    elseif type == "gq"
        order0 = γgq0_func(nf)
        order1 = γgq1_func(nf)
        order2 = γgq2_func(nf)
    elseif type == "gg"
        order0 = γgg0_func(nf)
        order1 = γgg1_func(nf)
        order2 = γgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end

#qq_x^2*log(x)-------------------------------------------------------------------------------------

function dγqq0_func(nf)

    if nf == 5
        value = 1.64335
    end

    return value
end

function dγqq1_func(nf)

    if nf == 5
        value = 6.30608
    end

    return value
end

function dγqq2_func(nf)

    if nf == 5
        value = 25.6848
    end

    return value
end

#qg_x^2*log(x)-------------------------------------------------------------------------------------

function dγqg0_func(nf)

    if nf == 5
        value = 0.661111
    end

    return value
end

function dγqg1_func(nf)

    if nf == 5
        value = -10.9396
    end

    return value
end

function dγqg2_func(nf)

    if nf == 5
        value = -25.7595
    end

    return value
end

#gq_x^2*log(x)-------------------------------------------------------------------------------------

function dγgq0_func(nf)

    if nf == 5
        value = 0.907407
    end

    return value
end

function dγgq1_func(nf)

    if nf == 5
        value = 15.1349
    end

    return value
end

function dγgq2_func(nf)

    if nf == 5
        value = -277.814
    end

    return value
end

#gg_x^2*log(x)-------------------------------------------------------------------------------------

function dγgg0_func(nf)

    if nf == 5
        value = 5.34254
    end

    return value
end

function dγgg1_func(nf)

    if nf == 5
        value = 31.7494
    end

    return value
end

function dγgg2_func(nf)

    if nf == 5
        value = -655.315
    end

    return value
end

function dγ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = dγqq0_func(nf)
        order1 = dγqq1_func(nf)
        order2 = dγqq2_func(nf)
    elseif type == "qg"
        order0 = dγqg0_func(nf)
        order1 = dγqg1_func(nf)
        order2 = dγqg2_func(nf)
    elseif type == "gq"
        order0 = dγgq0_func(nf)
        order1 = dγgq1_func(nf)
        order2 = dγgq2_func(nf)
    elseif type == "gg"
        order0 = dγgg0_func(nf)
        order1 = dγgg1_func(nf)
        order2 = dγgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end

#qq_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγqq0_func(nf)

    if nf == 5
        value = -0.541076
    end

    return value
end

function ddγqq1_func(nf)

    if nf == 5
        value = 2.33049
    end

    return value
end

function ddγqq2_func(nf)

    if nf == 5
        value = -35.2063
    end

    return value
end

#qg_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγqg0_func(nf)

    if nf == 5
        value = -0.435740
    end

    return value
end

function ddγqg1_func(nf)

    if nf == 5
        value = 10.7206
    end

    return value
end

function ddγqg2_func(nf)

    if nf == 5
        value = -39.8002
    end

    return value
end

#gq_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγgq0_func(nf)

    if nf == 5
        value = -1.02160
    end

    return value
end

function ddγgq1_func(nf)

    if nf == 5
        value = -0.135879
    end

    return value
end

function ddγgq2_func(nf)

    if nf == 5
        value = 777.820
    end

    return value
end

#gg_x^2*log(x)^2-----------------------------------------------------------------------------------

function ddγgg0_func(nf)

    if nf == 5
        value = -3.25458
    end

    return value
end

function ddγgg1_func(nf)

    if nf == 5
        value = 9.34699
    end

    return value
end

function ddγgg2_func(nf)

    if nf == 5
        value = 1472.63
    end

    return value
end

function ddγ_func(; type::String, α::Float64, order::Int64, nf::Int64)

    if     type == "qq"
        order0 = ddγqq0_func(nf)
        order1 = ddγqq1_func(nf)
        order2 = ddγqq2_func(nf)
    elseif type == "qg"
        order0 = ddγqg0_func(nf)
        order1 = ddγqg1_func(nf)
        order2 = ddγqg2_func(nf)
    elseif type == "gq"
        order0 = ddγgq0_func(nf)
        order1 = ddγgq1_func(nf)
        order2 = ddγgq2_func(nf)
    elseif type == "gg"
        order0 = ddγgg0_func(nf)
        order1 = ddγgg1_func(nf)
        order2 = ddγgg2_func(nf)
    else
        error("Type $type not supported.")
    end

    as = α/(4*π)

    if order < 0
        total = 0
    elseif order == 0
        total = as*order0
    elseif order == 1
        total = as*order0 + as^2*order1
    elseif order == 2
        total = as*order0 + as^2*order1 + as^3*order2
    end

    return total
end

#println(γqq_log_func(α=0.118, order=1, nf=5))
#println(γqg_log_func(α=0.118, order=1, nf=5))
#println(γgq_log_func(α=0.118, order=1, nf=5))
#println(γgg_log_func(α=0.118, order=1, nf=5))