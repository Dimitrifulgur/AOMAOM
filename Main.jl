using LinearAlgebra
using Combinatorics
using BenchmarkTools
using WignerSymbols
using PlotlyJS
using SparseArrays
using StaticArrays

include("Types.jl")
include("MathUtils.jl")
include("Integrals.jl")

function get_trigonal_prism(phi=0.0)
    ligands = [ 54.0  90.0+phi;
                54.0  210.0+phi;
                54.0  330.0+phi;
               126.0  90.0;
               126.0  210.0;
               126.0  330.0;]
    sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    return ligands, sigma
end

function get_octahedral(phi=0.0)
    ligands = [   0.0    0.0;
                 90.0    0.0;
                 90.0   90.0;
                 90.0  180.0;
                 90.0  270.0;
                180.0    0.0;]
    sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    return ligands, sigma
end

function get_trigonal_bipyramide(phi=0.0)
    ligands = [   0.0    0.0;
                 90.0   90.0;
                 90.0  210.0;
                 90.0  330.0;
                180.0    0.0;]
    sigma = [1000.0, 1000.0, 1000.0, 1000.0, 1000.0]
    return ligands, sigma
end

function get_tetrahedral(phi=0.0)
    theta = 110.0
    ligands = [   0.0    0.0;
                theta    0.0;
                theta  120.0;
                theta  240.0;]
    sigma = [1000.0, 1000.0, 1000.0, 1000.0]
    return ligands, sigma
end

function get_linear(phi=0.0)
    ligands = [   0.0    0.0+phi;
                180.0    0.0+phi;]
    sigma = [1000.0, 1000.0]
    return ligands, sigma
end

function main()
    nel = 7
    L = 2
    F0 = 194536.3
    F2 = 88404.7
    F4 = 55404.2
    
    # Выбираем геометрию (octahedral, tetrahedral, etc.)
    ligands, sigma = get_octahedral()
    
    # Считаем потенциал
    V = round.(AOM_system(ligands, sigma).V_AOM, digits = 5)
    V = AOMMatrixtoStandert(Matrix(V))
    V1 = RealtoYlm(ComplexF64.(V))
    
    Alldet = ConfigGen(nel, L)
    
    H = zeros(ComplexF64, size(Alldet,1), size(Alldet,1))
    H_Coulomb(Alldet, H, F0, F2, F4)
    Energ, EiVec = eigen(H)
    Energ = round.(Energ, digits = 20)
    EiVec = round.(EiVec, digits = 20)
    
    println("Кулоновские энергии (до внешних полей):")
    Energ_C = real.(Energ) .- minimum(real.(Energ))
    display(round.(Energ_C, digits=2))
    
    H = EiVec' * H * EiVec
    
    H_AOM(Alldet, EiVec, V1, H)
    H_SO(Alldet, EiVec, H, 120.0)
    
    En, Vn = eigen(H)
    Energ1 = real.(En) .- minimum(real.(En))
    Energ1 = round.(Energ1, digits=2)
    
    println("Расчитанные энергии:")
    display(Energ1)
end

# Запуск основной функции
main()