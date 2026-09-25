# ===============================
# From Integrals_Coulomb.jl
# ===============================
# <a11 b12 | e^2/rij |c21 d22>
# F0, F2, F4 - coulomb parameters
kr(x,y) = x == y ? 1.0 : 0.0
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
            # Coulomb operator conserves ML and MS
            if SlaterDeterminants[i].ml != SlaterDeterminants[j].ml || SlaterDeterminants[i].ms != SlaterDeterminants[j].ms
                continue
            end
            
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

# ===============================
# From Integrals_SO.jl
# ===============================
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

function H_SO(Alldet::Vector{SlaterDet}, coeff::Matrix{ComplexF64}, H::Matrix{ComplexF64}, λ::Float64)
    N = length(Alldet)
    H_det = zeros(ComplexF64, N, N)
    @inbounds for i in 1:N
        for j in i:N
            # Spin-orbit operator conserves M_J = M_L + M_S
            if (Alldet[i].ml + Alldet[i].ms) != (Alldet[j].ml + Alldet[j].ms)
                continue
            end
            
            val = OneElectronSODet(Alldet[i].det, Alldet[j].det, λ)
            H_det[i, j] = val
            H_det[j, i] = val'
        end
    end
    
    # H = H + coeff' * H_det * coeff
    mul!(H, coeff', H_det * coeff, 1.0, 1.0)
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

# ===============================
# From Integrals_AOM.jl
# ===============================
# такое чусвто, что здесь я оченб сильно увелчсиваю время расчета
# т.к. суммирую очень ненужные вклады в каждый элемент H[i,j]
# можно нахуярить еще матрицу между детерминантами лдя одноэлектронных интегралов
# чтобы каждый раз не смотреть перестановку
function H_AOM(Alldet::Vector{SlaterDet}, coeff::Matrix{ComplexF64}, V::Matrix{ComplexF64}, H::Matrix{ComplexF64})
    N = length(Alldet)
    H_det = zeros(ComplexF64, N, N)
    @inbounds for i in 1:N
        for j in i:N
            # Crystal field operators (like AOM) act only on spatial coordinates, conserving M_S
            if Alldet[i].ms != Alldet[j].ms
                continue
            end
            
            val = OneElectronAOMDet(Alldet[i].det, Alldet[j].det, V)
            H_det[i, j] = val
            H_det[j, i] = val'
        end
    end
    
    # H = H + coeff' * H_det * coeff
    mul!(H, coeff', H_det * coeff, 1.0, 1.0)
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

# ===============================
# From Integrals_Zeeman.jl
# ===============================
# такое чусвто, что здесь я оченб сильно увелчсиваю время расчета
# т.к. суммирую очень ненужные вклады в каждый элемент H[i,j]
# можно нахуярить еще матрицу между детерминантами лдя одноэлектронных интегралов
# чтобы каждый раз не смотреть перестановку
function H_Zeeman(Alldet::Vector{SlaterDet}, coeff::Matrix{ComplexF64}, H::Matrix{ComplexF64}, B::Vector{Float64})
    N = length(Alldet)
    H_det = zeros(ComplexF64, N, N)
    @inbounds for i in 1:N
        for j in i:N
            # Zeeman operator changes total M_J by at most 1
            if abs((Alldet[i].ml + Alldet[i].ms) - (Alldet[j].ml + Alldet[j].ms)) > 1.0
                continue
            end
            
            val = OneElectronZeemanDet(Alldet[i].det, Alldet[j].det, B)
            H_det[i, j] = val
            H_det[j, i] = val'
        end
    end
    
    mul!(H, coeff', H_det * coeff, 1.0, 1.0)
end


function OneElectronZeemanDet(det1::Vector{SpinOrbital}, det2::Vector{SpinOrbital}, B::Vector{Float64})
    # sort determinants and take sign of permutation
    sortdet = SortSlaterDet(det1, det2);
    det1_sorted = sortdet[1]; det2_sorted = sortdet[2]; sign_temp = sortdet[3]; ndif = sortdet[4];

    # take one electron integrals
    if ndif > 1
        return 0.0 + 0.0im
    elseif ndif == 1
        spinorb1 = det1_sorted[1]; spinorb2 = det2_sorted[1]; 
        integral = OneElectron_Zeeman(spinorb1, spinorb2, B);
        return sign_temp * integral
    else
        integral = 0.0 + 0.0im
        for i in axes(det1_sorted,1)
            integral += OneElectron_Zeeman(det1_sorted[i], det2_sorted[i], B);
        end
        return sign_temp * integral
    end
end

function OneElectron_Zeeman(SpinOrbital1::SpinOrbital, SpinOrbital2::SpinOrbital, B::Vector{Float64})
    Bx, By, Bz = B[1], B[2], B[3]
    # make hamiltonian in term basis
    l1 = SpinOrbital1.l; ml1 = SpinOrbital1.ml; s1 = SpinOrbital1.s; ms1 = SpinOrbital1.ms;
    l2 = SpinOrbital2.l; ml2 = SpinOrbital2.ml; s2 = SpinOrbital2.s; ms2 = SpinOrbital2.ms;
    kr1 = kr(l1, l2); kr2 = kr(ml1, ml2); kr3 = kr(s1,s2); kr4 = kr(ms1, ms2); kr5 = kr(ml1, (ml2+1));
    kr6 = kr((ms1+1), ms2);  kr7 = kr((ml1+1), ml2); kr8 = kr(ms1, (ms2+1));

    El = kr1*kr3*kr4*(kr5*(1 - im)/2*sqrt(l1*(l1+1) - ml1*ml2)*Bx + kr7*(1 + im)/2*sqrt(l1*(l1+1) - ml1*ml2)*Bx   + kr2*ml2*Bz);
    Es = kr1*kr2*kr3*(kr6*(1 - im)/2*sqrt(s1*(s1+1) - ms1*ms2)*Bx + kr8*(1 + im)/2*sqrt(s1*(s1+1) - ms1*ms2)*By   + kr2*ml2*Bz);
    Ezeem = El + Es;

    return Ezeem
end

