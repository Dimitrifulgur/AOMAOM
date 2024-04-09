# такое чусвто, что здесь я оченб сильно увелчсиваю время расчета
# т.к. суммирую очень ненужные вклады в каждый элемент H[i,j]
# можно нахуярить еще матрицу между детерминантами лдя одноэлектронных интегралов
# чтобы каждый раз не смотреть перестановку
function H_Zeeman(Alldet::Vector{SlaterDet}, coeff::Matrix{Float64}, H::Matrix{Basic})
    # I work in CSF basis
    # coeff - matrix of coefficients of CSF with respect to Slater Determinants 
    # ψ_i = coeff[:,i]
    # Можно брать не все коэффициенты!!!!
    # заранее найти только значимые и юзать их -> не буду считать 0 
    for i in axes(H,1)[begin:end]
        for j in axes(H,1)[begin:end]
            coeff1 = coeff[:,i];
            coeff2 = coeff[:,j];
            Energ =  CSF_Zeeman(Alldet, coeff1, coeff2);
             H[i,j] += Energ;
        end
    end
end

function CSF_Zeeman(Alldet::Vector{SlaterDet}, coeff1::Vector{Float64}, coeff2::Vector{Float64})
    E = 0.0;
    for i in eachindex(coeff1)
        for j in eachindex(coeff2)
            if (coeff1[i] == 0.0) || (coeff2[j] == 0.0)
                E += 0.0;
            else
                E += coeff1[i]*coeff2[j]*OneElectronZeemanDet(Alldet[i].det, Alldet[j].det);
            end
        end
    end
    return(E)
end


function OneElectronZeemanDet(Sdet1::SlaterDet, Sdet2::SlaterDet)
    # sort determinants and take sign of permutation
    sortdet = SortSlaterDet(Sdet1, Sdet2);
    det1 = sortdet[1]; det2 = sortdet[2]; sign_temp = sortdet[3]; ndif = sortdet[4];

    # take one electron integrals
    if ndif > 1
        return(0.0)
    elseif ndif == 1
        spinorb1 = det1[1]; spinorb2 = det2[1]; 
        integral = OneElectron_Zeeman(spinorb1, spinorb2);
        integral = sign_temp*integral;
        return(integral)
    else
        integral = 0.0
        for i in axes(det1,1)[begin:end]
            integral = integral + OneElectron_Zeeman(det1[i], det2[i]);
        end
        integral = sign_temp*integral;
        return(integral)
    end
end

function OneElectron_Zeeman(SpinOrbital1::SpinOrbital, SpinOrbital2::SpinOrbital)
    # make hamiltonian in term basis
    l1 = SpinOrbital1.l; ml1 = SpinOrbital1.ml; s1 = SpinOrbital1.s; ms1 = SpinOrbital1.ms;
    l2 = SpinOrbital2.l; ml2 = SpinOrbital2.ml; s2 = SpinOrbital2.s; ms2 = SpinOrbital2.ms;
    kr1 = kr(l1, l2); kr2 = kr(ml1, ml2); kr3 = kr(s1,s2); kr4 = kr(ms1, ms2); kr5 = kr(ml1, (ml2+1));
    kr6 = kr((ms1+1), ms2);  kr7 = kr((ml1+1), ml2); kr8 = kr(ms1, (ms2+1));

    El = kr1*kr3*kr4*(kr5*(1 - im)/2*sqrt(l1*(l1+1) - ml1*ml2)*Bx + kr7*(1 + im)/2*sqrt(l1*(l1+1) - ml1*ml2)*Bx   + kr2*ml2*Bz);
    Es = kr1*kr2*kr3*(kr6*(1 - im)/2*sqrt(s1*(s1+1) - ms1*ms2)*Bx + kr8*(1 + im)/2*sqrt(s1*(s1+1) - ms1*ms2)*By   + kr2*ml2*Bz);
    Ezeem = El + Es;

    return(Ezeem)

end