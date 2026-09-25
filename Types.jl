# ===============================
# From Type_SpinOrbital.jl
# ===============================
struct SpinOrbital
    l::Int
    ml::Int
    s::Float64
    ms::Float64
end

# ===============================
# From Type_SlaterDet.jl
# ===============================
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

# ===============================
# From Type_CSF.jl
# ===============================
struct CSF
    CoeffSlaterDet::Vector{Float64}
    l::Int
    s::Float64
    ml::Int
    ms::Float64
end

struct Term
    VecCSF::Vector{CSF}
    l::Int
    s::Float64
end

function Base.show(io::IO, t::Term)
    l = Int(t.l);
    s = t.s;
    TermSym = ["S", "P", "D", "F", "G", "H", "I", "K"]
    println(io, "$(Int(s*2+1))$(TermSym[l+1])\n")
end

 

# Check the purity of these CSFs with respect to spin and L
function DetermineTerm(E::Vector{Float64}, ψ::Matrix{Float64}, Alldet::Vector{SlaterDet})
    UniqueE = unique(E);
    Terms = Vector{Term}(undef, length(UniqueE))
    for i in eachindex(UniqueE)
        indexE = findall(x -> x == UniqueE[i], E);
        degenN = (indexE[end] - indexE[1] + 1);
        CSFvec = Vector{CSF}(undef, degenN);
        ml = zeros(Int, degenN);
        ms = zeros(Float64, degenN);

        for j in eachindex(indexE)
            NdetMax = findfirst(x -> x == maximum(abs.(ψ[:, indexE[j]])), abs.(ψ[:, indexE[j]]));
            ml[j] = Alldet[NdetMax].ml
            ms[j] = Alldet[NdetMax].ms
        end

        L = maximum(ml)
        S = maximum(ms)

        for j in eachindex(indexE)
            CSFvec[j] = CSF(ψ[:, indexE[j]], L, S, ml[j], ms[j])
        end

        Terms[i] = Term(CSFvec, L, S);
    end
    return Terms
end

function Ms(ψ::Vector{SlaterDet}, coeff::Vector{ComplexF64})
    ms = 0.0;
    for i in eachindex(coeff) 
        ms += abs2(coeff[i])*ψ[i].ms
    end
    return ms
end

function Ml(ψ::Vector{SlaterDet}, coeff::Vector{ComplexF64})
    ml = 0.0;
    for i in eachindex(coeff) 
        ml += abs2(coeff[i])*ψ[i].ml
    end
    return ml
end

