module JuliaInterface4CUTEstTest
  using Test, NLPModels, CUTEst
  include("./asserts/julia_interface_4_CUTEst.jl")

  using .JuliaInterface4CUTEst

  include("../algencan_interface.jl")

  problem_set::Vector{String} = CUTEst.select(
    min_var=20,max_var=50,objtype=[:quadratic,:sum_of_squares],custom_filter=x->x["constraints"]["ineq_both"] == 0
  )
  file = open("julia_interface_4_CUTEst_results.txt","w")

  for problem in problem_set
    
    @testset "Check algencan_interface behaviour over $problem problem" begin
      JuliaInterface4CUTEst.nlp = CUTEstModel(problem,lfirst=false,lvfirst=false)

      time_spent::Float64 = @elapsed begin
        x::Vector{Float64} = run_algencan()
      end
      
      f::Float64 = obj(JuliaInterface4CUTEst.nlp, x)

      println("========================RESULTS========================")
      println("f = $f")
      println("Time spent(seconds) = $time_spent")

      write(file,"$problem\t\t\t\t$time_spent\t\t\t\t$f\n")

      finalize(JuliaInterface4CUTEst.nlp)
    end
  end

  close(file)
end