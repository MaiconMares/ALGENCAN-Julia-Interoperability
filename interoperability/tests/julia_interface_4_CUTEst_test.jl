module JuliaInterface4CUTEstTest
  using Test, NLPModels, CUTEst
  include("./asserts/julia_interface_4_CUTEst.jl")

  using .JuliaInterface4CUTEst

  include("../algencan_interface.jl")

  problem_set::Vector{String} = CUTEst.select(max_var=10,custom_filter=x->x["constraints"]["ineq_both"] == 0)[1:10]

  for problem in problem_set
    
    @testset "Check algencan_interface behaviour over $problem problem" begin
      @time begin
        JuliaInterface4CUTEst.nlp = CUTEstModel(problem,lfirst=false,lvfirst=false)

        run_algencan()
        println("=======================TEST RESULTS=======================")
        print("Time spent: ")

        finalize(JuliaInterface4CUTEst.nlp)
      end
    end
  end
end