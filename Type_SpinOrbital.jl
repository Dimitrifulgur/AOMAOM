mutable struct SpinOrbital
    l::Number
    ml::Number
    s::Number
    ms::Number
    function SpinOrbital(l, ml, s, ms)
            new(l, ml, s, ms)
    end
end