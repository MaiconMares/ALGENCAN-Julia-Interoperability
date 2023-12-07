include("../tests/julia_interface_4_CUTEst_test.jl")
include("../tests/asserts/julia_interface_4_CUTEst.jl")
include("../problem_definition.jl")

using Documenter, CUTEst, .JuliaInterface4CUTEstTest, .JuliaInterface4CUTEst, .ProblemDefinition

makedocs(sitename="AlgencanInterface")