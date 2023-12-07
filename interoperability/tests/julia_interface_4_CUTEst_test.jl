"Defines a set of problems and run tests with them over ALGENCAN shared library."
module JuliaInterface4CUTEstTest
  using Test, NLPModels, CUTEst
  include("./asserts/julia_interface_4_CUTEst.jl")

  using .JuliaInterface4CUTEst

  include("../algencan_interface.jl")

  export run_tests

  """
    Executes tests over Julia interface with ALGENCAN from a defined set of 17 CUTEst problems.
  """
  function run_tests()::Nothing
    problem_set::Vector{String} = ["ANTWERP","HATFLDC","ORTHREGB","KSIP","BQPGABIM","OPTCNTRL","SANTALS","ERRINRSM","BQPGASIM","GOULDQP1","ACOPP14","ERRINROS","HATFLDGLS","TOINTQOR","3PK","METHANB8LS"]
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

  run_tests()
end