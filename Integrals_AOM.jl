# такое чусвто, что здесь я оченб сильно увелчсиваю время расчета
# т.к. суммирую очень ненужные вклады в каждый элемент H[i,j]
# можно нахуярить еще матрицу между детерминантами лдя одноэлектронных интегралов
# чтобы каждый раз не смотреть перестановку
function H_AOM(Alldet::Vector{SlaterDet}, coeff::Matrix{ComplexF64}, V::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    # I work in CSF basis
    # coeff - matrix of coefficients of CSF with respect to Slater Determinants 
    # ψ_i = coeff[:,i]
    # Можно брать не все коэффициенты!!!!
    # заранее найти только значимые и юзать их -> не буду считать 0 
    @inbounds for i in axes(H,1)
        @views coeff1 = coeff[:,i]
        for j in axes(H,1)[i:end]
            @views coeff2 = coeff[:,j]
            Energ = CSF_AOM(Alldet, coeff1, coeff2, V);
            if i != j  
                H[i,j] += Energ;
                H[j,i] += Energ';
            else
                H[i,j] +=  Energ;
            end     
        end
    end
end

function CSF_AOM(Alldet::Vector{SlaterDet}, coeff1::Vector{ComplexF64}, coeff2::Vector{ComplexF64}, V::Matrix{ComplexF64})
    E = 0.0;
    for i in eachindex(coeff1)
        for j in eachindex(coeff2)
            if (coeff1[i] == 0.0 + 0.0*im) || (coeff2[j] == 0.0 + 0.0*im)
                #E += 0.0;
                continue
                #println(i, "   ", j)
            else
                E += coeff1[i]*coeff2[j]*OneElectronAOMDet(Alldet[i].det, Alldet[j].det, V);
            end
        end
    end
    return(E)
end

function OneElectronAOMDet(Sdet1::Vector{SpinOrbital}, Sdet2::Vector{SpinOrbital}, V::Matrix{ComplexF64})
    # sort determinants and take sign of permutation
    sortdet = SortSlaterDet(Sdet1, Sdet2);
    det1 = sortdet[1]; det2 = sortdet[2]; sign_temp = sortdet[3]; ndif = sortdet[4];

    # take one electron integrals
    if ndif > 1
        return(0.0)
    elseif ndif == 1
        spinorb1 = det1[1]; spinorb2 = det2[1]; 
        integral = OneElectronAOM(spinorb1, spinorb2, V);
        integral = sign_temp*integral;
        return(integral)
    else
        integral = 0.0
        for i in axes(det1,1)[begin:end]
            integral = integral + OneElectronAOM(det1[i], det2[i], V);
        end
        integral = sign_temp*integral;
        return(integral)
    end
end

function OneElectronAOM(spinorbital1::SpinOrbital, spinorbital2::SpinOrbital, V::Matrix{ComplexF64})
    ms1 = spinorbital1.ms; ms2 = spinorbital2.ms;
    ml1 = spinorbital1.ml; ml2 = spinorbital2.ml;
    n = ml1 + 3; k = ml2  + 3; # ?????????
    delta = kr(ms1, ms2);
    # какой должен быть n,k ? V -> в базисе Y_lm
    if delta == 0.0
        return(0.0)
    else
        integral = delta*V[n, k];
    end
    return(integral)
end