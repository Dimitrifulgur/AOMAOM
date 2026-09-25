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

 

# проверить бы еще насколько эти CSF чистые по спину и L
function DetermineTerm(E::Vector{Float64}, ψ::Matrix{Float64}, Alldet::Vector{SlaterDet})
    UniqueE = unique(E);
    Terms = Array{Term}(undef, 1, size(UniqueE,1))
    for i in eachindex(UniqueE)
        indexE = findall(x -> x == UniqueE[i], E);
        degenN = (indexE[end] - indexE[1] + 1);
        CSFvec = Array{CSF}(undef, 1, degenN);
        ml = zeros(1, degenN);
        ms = zeros(1, degenN);

        for j in eachindex(indexE)
            NdetMax = findfirst(x -> x == maximum(abs.(ψ[:, indexE[j]])), abs.(ψ[:, indexE[j]])); # вот тут бы чекнуть по всем детерминантам какие ml и ms
            ml[j] = Alldet[NdetMax].ml
            ms[j] = Alldet[NdetMax].ms
        end

        L = maximum(ml)
        S = maximum(ms)

        for j in eachindex(indexE)
            CSFvec[j] = CSF(ψ[:, indexE[j]], L, S, ml[j], ms[j])
        end

        Terms[i] = Term(vec(CSFvec), L, S);
    end
    return(Terms)
end

function Ms(ψ::Vector{SlaterDet}, coeff::Vector{ComplexF64})
    ms = 0.0;
    for i in eachindex(coeff) 
        println(i, "      ", real(coeff[i]*coeff[i]'*ψ[i].ms))
        ms += coeff[i]*coeff[i]'*ψ[i].ms
    end
    return(ms)
end

function Ml(ψ::Vector{SlaterDet}, coeff::Vector{ComplexF64})
    ml = 0.0;
    for i in eachindex(coeff) 
        ml += coeff[i]*coeff[i]'*ψ[i].ml
    end
    return(ml)
end