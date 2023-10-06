module JuliaInterface4CUTEstTest
    using Test, NLPModels, CUTEst, LinearAlgebra
    include("./asserts/julia_interface_4_CUTEst.jl")

    using .JuliaInterface4CUTEst

    include("../algencan_interface.jl")
    JuliaInterface4CUTEst.nlp = CUTEstModel("HS27")

    @testset "Check algencan_interface behaviour over HS14 problem" begin
        @time begin
            #obj_func,x_vec = run_algencan()
            x = run_algencan()
            println("=======================TEST RESULTS=======================")
            print("Time spent: ")
        end
        #g = grad(nlp, x)
        #println(g)
        #gradient_norm = norm(g)
        #println("Objective function's gradient norm: $gradient_norm")
        
        finalize(JuliaInterface4CUTEst.nlp)
    end
end