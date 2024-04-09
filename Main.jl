include("./Included_Libs.jl")
#Base.ComplexF64(x::Basic) = convert(Complex{Float64}, x)
#conj(x::Basic) = real(x) ;


nel = 7;
L = 2;
# перегрузка для ComplexF64
F0 = 194536.3;
F2 = 88404.7;
F4 = 55404.2;
#=
#                    dz2       dxz        dyz       dx2-y2       dxy
sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]

global Energy = zeros(5,n+1);

for phi = 0:1:n
    ligands = [ 54.0  90.0+phi;
                54.0  210.0+phi;
                54.0  330.0+phi;
                126.0  90.0;
                126.0  210.0;
                126.0  330.0;]

    V = round.(AOM_system(ligands, sigma).V_AOM, digits = 15)
    V = AOMMatrixtoStandert(V)
    E,C = eigen(V)
    println(E)
    Energy[:, phi+1] = E;

end

phi = [i for i=0:1:n];
=#
#plot(phi, Energy')

#V = V - tr(V)/5*Matrix(I,5,5);

#= trigonal prism 74
phi = 0;
ligands = [ 54.0  90.0+phi;
    54.0  210.0+phi;
    54.0  330.0+phi;
    126.0  90.0;
    126.0  210.0;
    126.0  330.0;]
sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
=#

#= octahedral=   85.8
phi = 0;
ligands = [ 0.0  0.0;
            90.0  0.0;
            90.0  90.0;
            90.0  180.0;
            90.0  270.0;
            180.0  0.0;]
sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
=#

#=
#trigonal bipyramide  13.3 cm-1
phi = 0;
ligands = [ 0.0  0.0;
            90.0  90.0;
            90.0  210.0;
            90.0  330.0;
            180.0  0.0;]
sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
=#

#tetrahedral 1.5 cm-1
#=
phi = 0;
theta = 110.0
ligands = [ 0.0  0.0;
            theta  0.0;
            theta  120.0;
            theta  240.0;]
sigma = [1000.0, 1000.0, 1000.0, 1000.0]
=#

#linear 116.3
#=
phi = 0;
ligands = [ 0.0  0.0+phi;
            180.0  0.0+phi;]
sigma = [1000.0, 1000.0]
=#

V = round.(AOM_system(ligands, sigma).V_AOM, digits = 5)
V = AOMMatrixtoStandert(V)
V1 = RealtoYlm(ComplexF64.(V));
Alldet = ConfigGen(nel, L);
H = Matrix{ComplexF64}(undef, size(Alldet,1), size(Alldet,1));
H_Coulomb(Alldet, H, F0, F2, F4);
Energ, EiVec = eigen(H);
Energ = round.(Energ,digits = 20);
EiVec = round.(EiVec,digits = 20);
H = inv(EiVec)*H*EiVec;
#H = round.(H,digits=4)[1:28,1:28];
H_AOM(Alldet, EiVec, V1, H);
H_SO(Alldet, EiVec, H, 120.0);
En, Vn = eigen(H) ;
Energ1 = real.(En) .- minimum(real.(En));
Energ1 = round.(Energ1, digits=2);



#=
sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
global Energy = zeros(5,n+1);

for phi = 0:1:n
    ligands = [ 54.0  90.0+phi;
                54.0  210.0+phi;
                54.0  330.0+phi;
                126.0  90.0;
                126.0  210.0;
                126.0  330.0;]

    V = round.(AOM_system(ligands, sigma).V_AOM, digits = 15)
    V = AOMMatrixtoStandert(V)
    E,C = eigen(V)
    println(E)
    Energy[:, phi+1] = E;

end

phi = [i for i=0:1:n];E
=#