using Test
include("./asserts/hs11_problem_definition.jl")

using .HS11ProblemDefinition

include("../algencan_interface.jl")

@testset "Check algencan_interface behaviour over HS11 problem" begin
    @test 1 + 1 == 2
    @test 1 + 1 == 2
end;

run_algencan()