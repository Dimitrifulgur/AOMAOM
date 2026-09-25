struct SlaterDet
    det::Vector{SpinOrbital}
    ml::Int
    ms::Float64
    function SlaterDet(det::Vector{SpinOrbital})
        ml = 0
        ms = 0.0
        for i in det
            ml += i.ml
            ms += i.ms
        end
        new(det, ml, ms)
    end
end