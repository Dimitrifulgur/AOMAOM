mutable struct SlaterDet
    det::Vector{SpinOrbital}
    ml::Number
    ms::Number
    function SlaterDet(det::Vector{SpinOrbital})
        ml, ms = 0, 0;
        for i in det
            ml += i.ml;
            ms += i.ms;
        end
        new(det, ml, ms)
    end
end