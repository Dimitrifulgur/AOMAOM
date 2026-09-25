using LinearAlgebra
using Combinatorics
using BenchmarkTools
using WignerSymbols
using PlotlyJS
import Base
using SparseArrays
using StaticArrays
include("./Type_SpinOrbital.jl")
include("./Type_SlaterDet.jl")
include("./Type_CSF.jl")
include("./Functions_AOM.jl")

include("./Integrals_Coulomb.jl")
include("./Integrals_AOM.jl")
include("./Integrals_SO.jl")

include("./Other_SpinAlgebra.jl")
#include("./Integrals_Zeeman.jl")=#