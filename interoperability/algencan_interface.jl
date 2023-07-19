include("./problem_definition.jl")
using .ProblemDefinition
using Libdl
curr_dir = Vector{String}(["./"])
lib_path = Libdl.find_library("libalgencanma.so", curr_dir)
lib = Libdl.dlopen(lib_path, RTLD_NOW|RTLD_GLOBAL)

function run_algencan()
  evalf_ptr = @cfunction(ProblemDefinition.evalf!, Nothing, (Ref{Int32},Ptr{Float64},Ptr{Float64},Ref{Int32},Ref{MyDataPtr}))
  evalg_ptr = @cfunction(ProblemDefinition.evalg!, Nothing, (Ref{Int32},Ptr{Float64},Ptr{Float64},Ref{Int32},Ref{MyDataPtr}))
  evalc_ptr = @cfunction(ProblemDefinition.evalc!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Float64},Ref{Int32}, Ref{MyDataPtr}
  ))

  evalj_ptr = @cfunction(ProblemDefinition.evalj!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Int32},
    Ptr{Int32},Ptr{Int32},Ptr{Int32},Ref{Int32},
    Ptr{Int32},Ptr{Float64},Ref{Int32},Ref{MyDataPtr}
  ))
    
  evalhl_ptr = @cfunction(ProblemDefinition.evalhl!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Float64},Ref{Int32},
    Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},
    Ref{Int32},Ref{MyDataPtr}
  ))

  x,n,f,g,c,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin,inform,ind = problem_params()

  println("Calling function!")
  @ccall lib_path.__algencanma_MOD_init(
    evalf_ptr::Ptr{Cvoid}, evalg_ptr::Ptr{Cvoid},
    evalc_ptr::Ptr{Cvoid}, evalj_ptr::Ptr{Cvoid},
    evalhl_ptr::Ptr{Cvoid}, x::Ptr{Float64},
    n::Ref{Int32},f::Ref{Float64},g::Ptr{Float64},c::Ptr{Float64}, lind::Ptr{Int32}, lbnd::Ptr{Float64},
    uind::Ptr{Int32}, ubnd::Ptr{Float64},
    m::Ref{Int32}, p::Ref{Int32}, lambda::Ptr{Float64},
    jnnzmax::Ref{Int32}, hlnnzmax::Ref{Int32}, epsfeas::Ref{Float64},
    epscompl::Ref{Float64},epsopt::Ref{Float64}, rhoauto::Ref{Int32},
    rhoini::Ref{Float64},scale::Ref{Int32},extallowed::Ref{Int32},corrin::Ref{Int32},
    inform::Ref{Int32},ind::Ptr{Int32}
    )::Cvoid

  Libdl.dlclose(lib)
end