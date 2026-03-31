@everywhere begin

@inline function Σb_func(; as::T, Lb::T, Coeffs::AbstractVector{T}) where {T<:Real}
    @boundscheck @assert length(Coeffs) == 30

    # Destructure in same order (keeps existing meaning)
    (ΣbLL9NP,  ΣbNLL9NP,  ΣbNNLL9NP,
     ΣbLL8NP,  ΣbNLL8NP,  ΣbNNLL8NP,
     ΣbLL7NP,  ΣbNLL7NP,  ΣbNNLL7NP,
     ΣbLL6NP,  ΣbNLL6NP,  ΣbNNLL6NP,
     ΣbLL5NP,  ΣbNLL5NP,  ΣbNNLL5NP,
     ΣbLL4NP,  ΣbNLL4NP,  ΣbNNLL4NP,
     ΣbLL3NP,  ΣbNLL3NP,  ΣbNNLL3NP,
     ΣbLL2NP,  ΣbNLL2NP,  ΣbNNLL2NP,
     ΣbLL1NP,  ΣbNLL1NP,  ΣbNNLL1NP,
     ΣbLL0NP,  ΣbNLL0NP,  ΣbNNLL0NP) = Coeffs

    # Precompute powers once (Lb up to 9, as up to 9)
    Lb2 = Lb*Lb;   Lb3 = Lb2*Lb;   Lb4 = Lb3*Lb;   Lb5 = Lb4*Lb
    Lb6 = Lb5*Lb;  Lb7 = Lb6*Lb;   Lb8 = Lb7*Lb;   Lb9 = Lb8*Lb

    as2 = as*as;   as3 = as2*as;   as4 = as3*as;   as5 = as4*as
    as6 = as5*as;  as7 = as6*as;   as8 = as7*as;   as9 = as8*as

    # Terms (tower per order)
    t0 = ΣbLL0NP
    t1 = as  * (ΣbLL1NP*Lb + ΣbNLL1NP)
    t2 = as2 * (ΣbLL2NP*Lb2 + ΣbNLL2NP*Lb + ΣbNNLL2NP)
    t3 = as3 * (ΣbLL3NP*Lb3 + ΣbNLL3NP*Lb2 + ΣbNNLL3NP*Lb)
    t4 = as4 * (ΣbLL4NP*Lb4 + ΣbNLL4NP*Lb3 + ΣbNNLL4NP*Lb2)
    t5 = as5 * (ΣbLL5NP*Lb5 + ΣbNLL5NP*Lb4 + ΣbNNLL5NP*Lb3)
    t6 = as6 * (ΣbLL6NP*Lb6 + ΣbNLL6NP*Lb5 + ΣbNNLL6NP*Lb4)
    t7 = as7 * (ΣbLL7NP*Lb7 + ΣbNLL7NP*Lb6 + ΣbNNLL7NP*Lb5)
    t8 = as8 * (ΣbLL8NP*Lb8 + ΣbNLL8NP*Lb7 + ΣbNNLL8NP*Lb6)
    t9 = as9 * (ΣbLL9NP*Lb9 + ΣbNLL9NP*Lb8 + ΣbNNLL9NP*Lb7)

    return (((((((((t9 + t8) + t7) + t6) + t5) + t4) + t3) + t2) + t1) + t0)
end

#function Σb_func(; as::Float64, Lb::Float64, Coeffs::AbstractVector{<:Float64})
#
#    (ΣbLL9NP,  ΣbNLL9NP,  ΣbNNLL9NP,
#    ΣbLL8NP,  ΣbNLL8NP,  ΣbNNLL8NP,
#    ΣbLL7NP,  ΣbNLL7NP,  ΣbNNLL7NP,
#    ΣbLL6NP,  ΣbNLL6NP,  ΣbNNLL6NP,
#    ΣbLL5NP,  ΣbNLL5NP,  ΣbNNLL5NP,
#    ΣbLL4NP,  ΣbNLL4NP,  ΣbNNLL4NP,
#    ΣbLL3NP,  ΣbNLL3NP,  ΣbNNLL3NP,
#    ΣbLL2NP,  ΣbNLL2NP,  ΣbNNLL2NP,
#    ΣbLL1NP,  ΣbNLL1NP,  ΣbNNLL1NP,
#    ΣbLL0NP,  ΣbNLL0NP,  ΣbNNLL0NP) = Coeffs
#
#	ΣbNP = (
#		  as^0*(ΣbLL0NP                                      )
#		+ as^1*(ΣbLL1NP*Lb   + ΣbNLL1NP                      )
#		+ as^2*(ΣbLL2NP*Lb^2 + ΣbNLL2NP*Lb   + ΣbNNLL2NP     )
#		+ as^3*(ΣbLL3NP*Lb^3 + ΣbNLL3NP*Lb^2 + ΣbNNLL3NP*Lb  )
#		+ as^4*(ΣbLL4NP*Lb^4 + ΣbNLL4NP*Lb^3 + ΣbNNLL4NP*Lb^2)
#		+ as^5*(ΣbLL5NP*Lb^5 + ΣbNLL5NP*Lb^4 + ΣbNNLL5NP*Lb^3)
#		+ as^6*(ΣbLL6NP*Lb^6 + ΣbNLL6NP*Lb^5 + ΣbNNLL6NP*Lb^4)
#		+ as^7*(ΣbLL7NP*Lb^7 + ΣbNLL7NP*Lb^6 + ΣbNNLL7NP*Lb^5)
#		+ as^8*(ΣbLL8NP*Lb^8 + ΣbNLL8NP*Lb^7 + ΣbNNLL8NP*Lb^6)
#		+ as^9*(ΣbLL9NP*Lb^9 + ΣbNLL9NP*Lb^8 + ΣbNNLL9NP*Lb^7)
#	)
#	
#	return ΣbNP
#end

end