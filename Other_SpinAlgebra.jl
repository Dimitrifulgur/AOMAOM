function ckl(l::Number, k::Number, ma::Number, mc::Number)
    if abs(ma - mc) > k 
        return(0.)
    else
        ckl = (-1)^(ma)*(2*l+1)*wigner3j(l, k, l, 0, 0, 0)*wigner3j(l, k, l, -ma, ma - mc, mc)
        return(Float64(ckl))
    end
end


function SpinMatrices(S::Number)
    n = Int(2*S + 1)
    Sz =  zeros(ComplexF64, n, n); Sx =  zeros(ComplexF64, n, n); Sy =  zeros(ComplexF64, n, n); 
    Sp =  zeros(ComplexF64, n, n); Sm =  zeros(ComplexF64, n, n); S2 =  zeros(ComplexF64, n, n);
    for i in 1:n
        for j in 1:n
            mli = S - 1*(i-1)
            mlj = S - 1*(j-1)
            Sx[i, j] = (I[i, j-1] + I[i-1, j])/2*sqrt(S*(S+1)-mli*mlj)
            Sy[i, j] = (I[i, j-1] - I[i-1, j])/2/im*sqrt(S*(S+1)-mli*mlj)
            Sz[i, j] = I[i, j]*mli
            Sp[i, j] = I[i, j-1]*sqrt(S*(S+1)-mli*mlj)
            Sm[i, j] = I[i-1, j]*sqrt(S*(S+1)-mli*mlj)
            S2[i, j] = I[i, j]*sqrt(S*(S+1))
        end
    end
    return(Sx, Sy, Sz, Sp, Sm, S2)
end

function RealtoYlm(V::Matrix{ComplexF64})
    # R is the transormation matrix from (z2 xz yz x2-y2 xy) to (-2 -1  0  1  2)
    R = [0.0        0.0            0.0       1/sqrt(2)     -im/sqrt(2);
         0.0     1/sqrt(2)     -im/sqrt(2)     0.0             0.0    ;
         1.0        0.0           0.0         0.0             0.0    ;
         0.0     -1/sqrt(2)   -im/sqrt(2)     0.0             0.0    ;
         0.0        0.0           0.0       1/sqrt(2)      im/sqrt(2) ;];

    R1 = [0.0           0.0           1.0       0.0             0.0;
          0.0         1/sqrt(2)       0.0      -1/sqrt(2)       0.0    ;
          0.0         im/sqrt(2)       0.0      im/sqrt(2)        0.0    ;
          1/sqrt(2)     0.0           0.0        0.0          1/sqrt(2)  ;
          im/sqrt(2)      0.0          0.0        0.0          -im/sqrt(2) ;];
    Vnew = R*V*R1;
    return(Vnew)
end

function ConfigGen(n::Number, l::Number)
    ndet = Int(factorial((2*l + 1)*2)/factorial(n)/factorial((2*l+1)*2 - n));
    det = Vector{SlaterDet}(undef, ndet)
    ml = [i for i in -l:l];
    ms = [i for i in -0.5:1:0.5];
    dset1 = [SpinOrbital(l, i, 1/2, ms[1]) for i in ml];
    dset2 = [SpinOrbital(l, i, 1/2, ms[2]) for i in ml];
    dset = [dset1; dset2];
    det_temp = collect(combinations(dset, n));

    for i in det_temp
        det = [SlaterDet(i) for i in det_temp]
    end
    return(det)
end

# p is a degree in 10^(-p)
# p = presicion
#=function SimplifyRound(t::Basic, p::Int)
    symbolt = Symbolics.get_variables(t);
    tnew = Basic(0.0);
    for i in eachindex(symbolt)
        tnew += Basic(round(Symbolics.coeff(t, symbolt[i]), digits=p)*symbolt[i]);
    end
    return(tnew)
end

function Base.:(==)(a::Basic, b::Basic)
    symbolA = Symbolics.get_variables(a);
    symbolB = Symbolics.get_variables(b);
    keyBool = true;
    if size(symbolA, 1) == size(symbolB, 1)
        for i in eachindex(symbolA)
            keyBool = keyBool && (Symbolics.coeff(a, symbolA[i]) == Symbolics.coeff(b, symbolB[i]))
        end
        return(keyBool)
    else
        return(false)
    end
end

function Base.:(!=)(a::Basic, b::Basic)
    return(!(a==b))
end
# wigner3j(l, k, l, 0, 0, 0)
=#