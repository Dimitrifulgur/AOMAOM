# z2 yz xz xy x2-y2
#########################
#############  σ function  ##################
#10.1021/acs.inorgchem.5b01706
Fσz2(θ, φ) = 1/4*(1 + 3*cos(2*θ));
Fσyz(θ, φ) = sqrt(3)/2*sin(φ)*sin(2*θ);
Fσxz(θ, φ) = sqrt(3)/2*cos(φ)*sin(2*θ);
Fσxy(θ, φ) = sqrt(3)/4*sin(2*φ)*(1.0 - cos(2*θ));
Fσx2y2(θ, φ) = sqrt(3)/4*cos(2*φ)*(1.0 - cos(2*θ));

#############  πs function  ##################
Fπsz2(θ, φ, ψ) = sqrt(3)/2*sin(2*θ)*sin(ψ);
Fπsyz(θ, φ, ψ) = cos(φ)*cos(θ)*cos(ψ) - sin(φ)*cos(2*θ)*sin(ψ);
Fπsxz(θ, φ, ψ) = -sin(φ)*cos(θ)*cos(ψ) - cos(φ)*cos(2*θ)*sin(ψ);
Fπsxy(θ, φ, ψ) = cos(2*φ)*sin(θ)*cos(ψ) - 1/2*sin(2*φ)*sin(2*θ)*sin(ψ);
Fπsx2y2(θ, φ, ψ) = -sin(2*φ)*sin(θ)*cos(ψ) - 1/2*cos(2*φ)*sin(2*θ)*sin(ψ);

#############  πc function  ##################
Fπcz2(θ, φ, ψ) = -sqrt(3)/2*sin(2*θ)*cos(ψ);
Fπcyz(θ, φ, ψ) = cos(φ)*cos(θ)*sin(ψ) + sin(φ)*cos(2*θ)*cos(ψ);
Fπcxz(θ, φ, ψ) = -sin(φ)*cos(θ)*sin(ψ) + cos(φ)*cos(2*θ)*cos(ψ);
Fπcxy(θ, φ, ψ) = cos(2*φ)*sin(θ)*sin(ψ) + 1/2*sin(2*φ)*sin(2*θ)*cos(ψ);
Fπcx2y2(θ, φ, ψ) = -sin(2*φ)*sin(θ)*sin(ψ) + 1/2*cos(2*φ)*sin(2*θ)*cos(ψ);

Fσπcz2(θ, φ, ψ) = -sqrt(3)*(3*cos(θ)^2-1)*sin(2*θ);
Fσπcyz(θ, φ, ψ) =  sqrt(3)/4*sin(2*θ)*cos(2*θ);
Fσπcxz(θ, φ, ψ) =  0.0;
Fσπcxy(θ, φ, ψ) =  sqrt(3)/8*(1-cos(2*θ)*sin(2*θ));
Fσπcx2y2(θ, φ, ψ) = 0.0;

function SigmaAOM(θ::Float64, φ::Float64)
    θ = θ/180*π;  φ = φ/180*π;
    
    V11 = Fσz2(θ, φ)*Fσz2(θ, φ); V12 = Fσz2(θ, φ)*Fσyz(θ, φ); V13 = Fσz2(θ, φ)*Fσxz(θ, φ);  V14 = Fσz2(θ, φ)*Fσxy(θ, φ);   V15 = Fσz2(θ, φ)*Fσx2y2(θ, φ);
    V22 = Fσyz(θ, φ)*Fσyz(θ, φ); V23 = Fσyz(θ, φ)*Fσxz(θ, φ); V24 = Fσyz(θ, φ)*Fσxy(θ, φ);  V25 = Fσyz(θ, φ)*Fσx2y2(θ, φ);
    V33 = Fσxz(θ, φ)*Fσxz(θ, φ); V34 = Fσxz(θ, φ)*Fσxy(θ, φ); V35 = Fσxz(θ, φ)*Fσx2y2(θ, φ);
    V44 = Fσxy(θ, φ)*Fσxy(θ, φ); V45 = Fσxy(θ, φ)*Fσx2y2(θ, φ);
    V55 = Fσx2y2(θ, φ)*Fσx2y2(θ, φ);

    Vσ_matrix = [V11 V12 V13 V14 V15; 
                 V12 V22 V23 V24 V25; 
                 V13 V23 V33 V34 V35;
                 V14 V24 V34 V44 V45;
                 V15 V25 V35 V45 V55];

    return(Vσ_matrix)
end

function PiSineAOM(θ::Float64, φ::Float64, ψ::Float64)
    θ = θ/180*π;  φ = φ/180*π; ψ = ψ/180*π;

    V11 = Fπsz2(θ, φ, ψ)*Fπsz2(θ, φ, ψ); V12 = Fπsz2(θ, φ, ψ)*Fπsyz(θ, φ, ψ); V13 = Fπsz2(θ, φ, ψ)*Fπsxz(θ, φ, ψ);   V14 = Fπsz2(θ, φ, ψ)*Fπsxy(θ, φ, ψ);   V15 = Fπsz2(θ, φ, ψ)*Fπsx2y2(θ, φ, ψ);
    V22 = Fπsyz(θ, φ, ψ)*Fπsyz(θ, φ, ψ); V23 = Fπsyz(θ, φ, ψ)*Fπsxz(θ, φ, ψ); V24 = Fπsyz(θ, φ, ψ)*Fπsxy(θ, φ, ψ);   V25 = Fπsyz(θ, φ, ψ)*Fπsx2y2(θ, φ, ψ); V55 = Fπsx2y2(θ, φ, ψ)*Fπsx2y2(θ, φ, ψ);
    V33 = Fπsxz(θ, φ, ψ)*Fπsxz(θ, φ, ψ); V34 = Fπsxz(θ, φ, ψ)*Fπsxy(θ, φ, ψ); V35 = Fπsxz(θ, φ, ψ)*Fπsx2y2(θ, φ, ψ); V44 = Fπsxy(θ, φ, ψ)*Fπsxy(θ, φ, ψ);   V45 = Fπsxy(θ, φ, ψ)*Fπsx2y2(θ, φ, ψ);

    Vπs_matrix = [V11 V12 V13 V14 V15; 
                 V12 V22 V23 V24 V25; 
                 V13 V23 V33 V34 V35;
                 V14 V24 V34 V44 V45;
                 V15 V25 V35 V45 V55];

    return(Vπs_matrix)
end

function PiCosineAOM(θ::Float64, φ::Float64, ψ::Float64)
    θ = θ/180*π;  φ = φ/180*π; ψ = ψ/180*π;

    V11 = Fπcz2(θ, φ, ψ)*Fπcz2(θ, φ, ψ); V12 = Fπcz2(θ, φ, ψ)*Fπcyz(θ, φ, ψ); V13 = Fπcz2(θ, φ, ψ)*Fπcxz(θ, φ, ψ);   V14 = Fπcz2(θ, φ, ψ)*Fπcxy(θ, φ, ψ);   V15 = Fπcz2(θ, φ, ψ)*Fπcx2y2(θ, φ, ψ);
    V22 = Fπcyz(θ, φ, ψ)*Fπcyz(θ, φ, ψ); V23 = Fπcyz(θ, φ, ψ)*Fπcxz(θ, φ, ψ); V24 = Fπcyz(θ, φ, ψ)*Fπcxy(θ, φ, ψ);   V25 = Fπcyz(θ, φ, ψ)*Fπcx2y2(θ, φ, ψ); V55 = Fπcx2y2(θ, φ, ψ)*Fπcx2y2(θ, φ, ψ);
    V33 = Fπcxz(θ, φ, ψ)*Fπcxz(θ, φ, ψ); V34 = Fπcxz(θ, φ, ψ)*Fπcxy(θ, φ, ψ); V35 = Fπcxz(θ, φ, ψ)*Fπcx2y2(θ, φ, ψ); V44 = Fπcxy(θ, φ, ψ)*Fπcxy(θ, φ, ψ);   V45 = Fπcxy(θ, φ, ψ)*Fπcx2y2(θ, φ, ψ);

    Vπs_matrix = [V11 V12 V13 V14 V15; 
                 V12 V22 V23 V24 V25; 
                 V13 V23 V33 V34 V35;
                 V14 V24 V34 V44 V45;
                 V15 V25 V35 V45 V55];

    return(Vπs_matrix)
end

function PiCosineSigmaAOM(θ::Float64, φ::Float64, ψ::Float64)
    θ = θ/180*π;  φ = φ/180*π; ψ = ψ/180*π;
    V11 = Fσπcz2(θ, φ, ψ)*Fσπcz2(θ, φ, ψ); V12 = Fσπcz2(θ, φ, ψ)*Fσπcyz(θ, φ, ψ); V13 = Fσπcz2(θ, φ, ψ)*Fσπcxz(θ, φ, ψ);   V14 = Fσπcz2(θ, φ, ψ)*Fσπcxy(θ, φ, ψ);   V15 = Fσπcz2(θ, φ, ψ)*Fσπcx2y2(θ, φ, ψ);
    V22 = Fσπcyz(θ, φ, ψ)*Fσπcyz(θ, φ, ψ); V23 = Fσπcyz(θ, φ, ψ)*Fσπcxz(θ, φ, ψ); V24 = Fσπcyz(θ, φ, ψ)*Fσπcxy(θ, φ, ψ);   V25 = Fσπcyz(θ, φ, ψ)*Fσπcx2y2(θ, φ, ψ); V55 = Fσπcx2y2(θ, φ, ψ)*Fσπcx2y2(θ, φ, ψ);
    V33 = Fσπcxz(θ, φ, ψ)*Fσπcxz(θ, φ, ψ); V34 = Fσπcxz(θ, φ, ψ)*Fσπcxy(θ, φ, ψ); V35 = Fσπcxz(θ, φ, ψ)*Fσπcx2y2(θ, φ, ψ); V44 = Fσπcxy(θ, φ, ψ)*Fσπcxy(θ, φ, ψ);   V45 = Fσπcxy(θ, φ, ψ)*Fσπcx2y2(θ, φ, ψ);

    Vπσs_matrix = [V11 V12 V13 V14 V15; 
                 V12 V22 V23 V24 V25; 
                 V13 V23 V33 V34 V35;
                 V14 V24 V34 V44 V45;
                 V15 V25 V35 V45 V55];

    return(Vπσs_matrix)
end


mutable struct AOM_system
    ligands::Matrix{Float64}     # [θ1, φ1, ψ1; θ2, φ2, ψ2; θ3, φ3, ψ3; ..... ]
    eσ::Vector{Float64}
    V_AOM::Matrix{Float64}
    function AOM_system(ligands::Matrix{Float64}, sigma::Vector{Float64})
        V_AOM = zeros(Float64, 5, 5);
        for i in 1:size(ligands,1)
            θ, φ  = ligands[i, :];
            V_AOM += sigma[i]*SigmaAOM(θ, φ);
        end
        new(ligands, sigma, V_AOM)
    end
end


###################################################################################################################################################
function AOMMatrixtoStandert(V::Matrix{Float64})
    Vnew = Matrix{Float64}(undef,5,5);
    Vnew[1,1] = V[1,1]; Vnew[2,1] = Vnew[1,2] = V[3,1]; Vnew[3,1] = Vnew[1,3] = V[2,1]; Vnew[4,1] = Vnew[1,4] = V[5,1] ; Vnew[5,1] = Vnew[1,5] = V[4,1];
    Vnew[2,2] = V[3,3]; Vnew[3,2] = Vnew[2,3] = V[2,3]; Vnew[4,2] = Vnew[2,4] = V[3,5]; Vnew[5,2] = Vnew[2,5] = V[4,3] ;
    Vnew[3,3] = V[2,2]; Vnew[4,3] = Vnew[3,4] = V[5,2]; Vnew[5,3] = Vnew[3,5] = V[4,2];
    Vnew[4,4] = V[5,5]; Vnew[4,5] = Vnew[5,4] = V[5,4]; 
    Vnew[5,5] = V[4,4]; 

    return(Vnew)
end
 #=
#        dz2      dxz     dyz     dx2-y2     dxy
# dz2
# dxz
# dyz     Vnew            
# dx2-y2
# dxy
##############################################
#         z2    yz    xz   xy    x2-y2  
# z2
# yz     V       
# xz             
# xy
# dx2y2

###################################################################################################################################################
=#