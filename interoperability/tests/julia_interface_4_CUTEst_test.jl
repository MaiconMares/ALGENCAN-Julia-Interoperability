module JuliaInterface4CUTEstTest
    using Test, NLPModels, CUTEst, LinearAlgebra
    include("./asserts/julia_interface_4_CUTEst.jl")

    using .JuliaInterface4CUTEst

    include("../algencan_interface.jl")

    problem_set::Vector{String} = CUTEst.select(max_var=10,custom_filter=x->x["constraints"]["ineq_both"] == 0)[1:3]
    #problem_set::Vector{String} = ["POLAK4","EXPFITA"]

    for problem in problem_set
        JuliaInterface4CUTEst.nlp = CUTEstModel(problem,lfirst=false,lvfirst=false)

        @testset "Check algencan_interface behaviour over $problem problem" begin
            @time begin
                run_algencan()
                println("=======================TEST RESULTS=======================")
                print("Time spent: ")
                finalize(JuliaInterface4CUTEst.nlp)
            end
        end
    end
end