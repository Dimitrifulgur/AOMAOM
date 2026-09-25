# <a11 b12 | e^2/rij |c21 d22>
# F0, F2, F4 - coulomb parameters
kr(x,y) = ==(x,y)*1.0;
function TwoElectronCoulomb(spinorb11::SpinOrbital, spinorb12::SpinOrbital, spinorb21::SpinOrbital, spinorb22::SpinOrbital, F0::Float64, F2::Float64, F4::Float64)
    delta1 = kr(spinorb11.ms, spinorb21.ms);
    delta2 = kr(spinorb12.ms, spinorb22.ms);
    delta3 = kr((spinorb11.ml + spinorb12.ml),(spinorb21.ml + spinorb22.ml));
    l = spinorb11.l;
    ml11 = spinorb11.ml;  ml12 = spinorb12.ml;
    ml21 = spinorb21.ml;  ml22 = spinorb22.ml;  

    #delta1 * delta2 * delta3 * (F0*ckl(l, 0, ml11, ml21)*ckl(l, 0, ml22, ml12) + F2*ckl(l, 2, ml11, ml21)*ckl(l, 2, ml22, ml12) + F4*ckl(l, 4, ml11, ml21)*ckl(l, 4, ml22, ml12))
    if (delta1 == 0.0) || (delta2 == 0.0) ||(delta3 == 0.0)
        return(0.0)
    else
      #  integral = delta1 * delta2 * delta3 * (F0*ckl(l, 0, ml11, ml21)*ckl(l, 0, ml22, ml12) + 49*F2*ckl(l, 2, ml11, ml21)*ckl(l, 2, ml22, ml12) + 441*F4*ckl(l, 4, ml11, ml21)*ckl(l, 4, ml22, ml12))
      integral = delta1 * delta2 * delta3 * (F0*ckl(l, 0, ml11, ml21)*ckl(l, 0, ml22, ml12) + F2*ckl(l, 2, ml11, ml21)*ckl(l, 2, ml22, ml12) + F4*ckl(l, 4, ml11, ml21)*ckl(l, 4, ml22, ml12))  
      return(integral)
    end
end

# function to turn Slater determinants into standart order
function SortSlaterDet(det11::Vector{SpinOrbital}, det22::Vector{SpinOrbital})
    det1 = copy(det11); det2 = copy(det22);
    accord1in2 = [x in det2 for x = det1];
    accord2in1 = [x in det1 for x = det2];
    n_match = sum([x in det1 for x = det2]);
    n_dismatch = length(accord1in2) - n_match;
    sign_perm1 = 1; # sign if permutation
    sign_perm2 = 1;

    # different spin-orbitals in beginning of the vector
    m = 1;
    for i in axes(accord1in2,1)[begin:end]
        if (accord1in2[i] == 0)
            det1[i], det1[m] = det1[m], det1[i];
            if i == m
                sign_perm1 = sign_perm1;
                m += 1;
            else
                sign_perm1 = -sign_perm1;
                m += 1;
            end
        end
    end

    j = 1;
    for i in axes(accord2in1,1)[begin:end]
        if (accord2in1[i] == 0)
            det2[i], det2[j] = det2[j], det2[i];
            if i == j
                sign_perm2 = sign_perm2;
                j += 1;
            else
                sign_perm2 = -sign_perm2;
                j += 1;
            end
        end
    end

    # sort matching part of the vector
    # take det1 and find permuation sign for  matching part in det2
    for i in (n_dismatch+1):length(det1)
        for k in i:length(det1)
            if (det1[i] == det2[k]) & (k != i)
                det2[i], det2[k] = det2[k], det2[i];
                sign_perm2 = -sign_perm2;
            end
        end
    end

    return(det1, det2, sign_perm1*sign_perm2, n_dismatch)    
end

# take Coulomb integrals between two Slater integrals
function CoulombIntegralDet(det11::Vector{SpinOrbital}, det22::Vector{SpinOrbital}, F0::Float64, F2::Float64, F4::Float64)
    # sort determinants and take sigh of permutation
    sortdet = SortSlaterDet(det11, det22);
    det1 = sortdet[1]; det2 = sortdet[2]; sign_temp = sortdet[3]; ndif = sortdet[4];

    if ndif > 2
        return(0.)

    elseif ndif == 2
        # take different orbitals in determinants
        spinorb11 = det1[1]; spinorb12 = det1[2]; 
        spinorb21 = det2[1]; spinorb22 = det2[2];

        coulomb = TwoElectronCoulomb(spinorb11, spinorb12, spinorb21, spinorb22, F0, F2, F4);
        exchange = TwoElectronCoulomb(spinorb11, spinorb12, spinorb22, spinorb21, F0, F2, F4);
        integral = sign_temp*(coulomb - exchange);
        return(integral)
    elseif ndif == 1
        # take different orbitals in determinants
        spinorb11 = det1[1]; spinorb21 = det2[1]; 
        integral = 0.;
        for i in axes(det11,1)[2:end]
            spinorbit_i = det1[i]
            coulomb = TwoElectronCoulomb(spinorb11, spinorbit_i, spinorb21, spinorbit_i, F0, F2, F4);
            exchange = TwoElectronCoulomb(spinorb11, spinorbit_i, spinorbit_i, spinorb21, F0, F2, F4);
            integral = integral + coulomb - exchange;
        end
        integral = sign_temp*integral;
        return(integral)

    elseif ndif == 0
        integral = 0.;
        for i in axes(det11,1)[1:end]
            for j in axes(det11,1)[i:end]
                spinorb11 = det1[i]; spinorb12 = det1[j]; 
                coulomb = TwoElectronCoulomb(spinorb12, spinorb11, spinorb12, spinorb11, F0, F2, F4);
                exchange = TwoElectronCoulomb(spinorb12, spinorb11, spinorb11, spinorb12, F0, F2, F4);
                integral = integral + coulomb - exchange;
            end
        end
        integral = sign_temp*integral;
        return(integral)
    end
end


function H_Coulomb(SlaterDeterminants::Vector{SlaterDet}, H::Matrix{ComplexF64}, F0::Float64, F2::Float64, F4::Float64)
    Alldet = [SlaterDeterminants[i].det for i in eachindex(SlaterDeterminants)]
    @inbounds for i in axes(Alldet,1)[begin:end]
        for j in axes(Alldet,1)[i:end]
            # опасный момент
            if i == j
                H[i,j] = CoulombIntegralDet(Alldet[i], Alldet[j], F0, F2, F4);
            else
                val = CoulombIntegralDet(Alldet[i], Alldet[j], F0, F2, F4);
                H[i,j] = val;
                H[j,i] = val'
            end
        end
    end
    return(1)
end