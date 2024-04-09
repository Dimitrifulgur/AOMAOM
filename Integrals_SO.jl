#=function HTerm_SO(H::Matrix{Num}, Terms::Matrix{Term})
    CSF_full = Terms[1].VecCSF[begin:end]

    for i in eachindex(Terms)
        if i != 1
            append!(CSF_full, Terms[i].VecCSF[begin:end])
        end
    end

    for i in eachindex(CSF_full)
        for j in eachindex(CSF_full)
            H[i, j] += OneElectronSODet(CSF_full[i], CSF_full[j])
        end
    end
end
=#

# такое чусвто, что здесь я оченб сильно увелчсиваю время расчета
# т.к. суммирую очень ненужные вклады в каждый элемент H[i,j]
# можно нахуярить еще матрицу между детерминантами лдя одноэлектронных интегралов
# чтобы каждый раз не смотреть перестановку
function H_SO(Alldet::Vector{SlaterDet}, coeff::Matrix{ComplexF64}, H::Matrix{ComplexF64}, λ::Float64)
    # I work in CSF basis
    # coeff - matrix of coefficients of CSF with respect to Slater Determinants 
    # ψ_i = coeff[:,i]
    # Можно брать не все коэффициенты!!!!
    # заранее найти только значимые и юзать их -> не буду считать 0 
    for i in axes(H,1)[begin:end]
        for j in axes(H,1)[begin:end]
            coeff1 = coeff[:,i];
            coeff2 = coeff[:,j];
            Energ = CSF_SO(Alldet, coeff1, coeff2, λ);
            H[i,j] +=  Energ ;
        end
    end
end

function CSF_SO(Alldet::Vector{SlaterDet}, coeff1::Vector{ComplexF64}, coeff2::Vector{ComplexF64}, λ::Float64)
    E = 0.0;
    for i in eachindex(coeff1)
        for j in eachindex(coeff2)
            if (coeff1[i] == 0.0) || (coeff2[j] == 0.0)
                E += 0.0;
            else
                E += coeff1[i]*coeff2[j]*OneElectronSODet(Alldet[i].det, Alldet[j].det, λ);
            end
        end
    end
    return(E)
end


function OneElectronSODet(Sdet1::Vector{SpinOrbital}, Sdet2::Vector{SpinOrbital}, λ::Float64)
    # sort determinants and take sign of permutation
    sortdet = SortSlaterDet(Sdet1, Sdet2);
    det1 = sortdet[1]; det2 = sortdet[2]; sign_temp = sortdet[3]; ndif = sortdet[4];


    # take one electron integrals
    if ndif > 1
        return(0.0)
    elseif ndif == 1
        spinorb1 = det1[1]; spinorb2 = det2[1]; 
        integral = OneElectron_SO(spinorb1, spinorb2, λ);
        integral = sign_temp*integral;
        return(integral)
    else
        integral = 0.0
        for i in axes(det1,1)[begin:end]
            integral = integral + OneElectron_SO(det1[i], det2[i], λ);
        end
        integral = sign_temp*integral;
        return(integral)
    end
end

function OneElectron_SO(SpinOrbital1::SpinOrbital, SpinOrbital2::SpinOrbital, λ::Float64)
    # make hamiltonian in term basis
    l1 = SpinOrbital1.l; ml1 = SpinOrbital1.ml; s1 = SpinOrbital1.s; ms1 = SpinOrbital1.ms;
    l2 = SpinOrbital2.l; ml2 = SpinOrbital2.ml; s2 = SpinOrbital2.s; ms2 = SpinOrbital2.ms;
    kr1 = kr(l1, l2); kr2 = kr(ml1, ml2); kr3 = kr(s1,s2); kr4 = kr(ms1, ms2); kr5 = kr(ml1, (ml2+1));
    kr6 = kr((ms1+1), ms2);  kr7 = kr((ml1+1), ml2); kr8 = kr(ms1, (ms2+1));

    ESO = λ*kr1*kr3*(ml1*ms1*kr2*kr4 + 1/2*sqrt(l1*(l1+1) - ml1*ml2)*sqrt(s1*(s1+1) - ms1*ms2)*(kr5*kr6 + kr7*kr8));
    return(ESO)
end