include("./problem_definition.jl")
using .ProblemDefinition
using Libdl
curr_dir = Vector{String}(["./"])
lib_path = Libdl.find_library("libalgencanma.so", curr_dir)
lib = Libdl.dlopen(lib_path)

function run_algencan()
  evalf_ptr = @cfunction(ProblemDefinition.evalf!, Nothing, (Ref{Int64},Ptr{Float64},Ptr{Float64},Ptr{Int64}))
  evalg_ptr = @cfunction(ProblemDefinition.evalg!, Nothing, (Ref{Int64},Ptr{Float64},Ptr{Float64},Ref{Int64}))
  evalc_ptr = @cfunction(ProblemDefinition.evalc!, Nothing, (
    Ref{Int64},Ptr{Float64},Ref{Int64},Ref{Int64},Ptr{Float64},Ptr{Int64}
  ))

  evalj_ptr = @cfunction(ProblemDefinition.evalj!, Nothing, (
    Ref{Int64},Ptr{Float64},Ref{Int64},Ref{Int64},Ref{Vector{Int32}},
    Ref{Vector{Int32}},Ref{Vector{Int64}},Ref{Vector{Int64}},Ref{Int64},
    Ref{Vector{Int64}},Ptr{Float64},Ref{Int64},Ref{MyDataPtr}
  ))

  evalhl_ptr = @cfunction(ProblemDefinition.evalhl!, Nothing, (
    Ref{Int64},Ptr{Float64},Ref{Int64},Ref{Int64},Ptr{Float64},Ref{Int64},
    Ref{Int32},Ref{Int64},Ref{Vector{Int64}},Ref{Vector{Int64}},Ptr{Float64},
    Ref{Int64},Ref{Int64},Ref{MyDataPtr}
  ))

  x,n,f,g,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin,inform = problem_params()


  @ccall lib_path.__algencanma_MOD_init(
    evalf_ptr::Ptr{Cvoid}, evalg_ptr::Ptr{Cvoid},
    evalc_ptr::Ptr{Cvoid}, evalj_ptr::Ptr{Cvoid},
    evalhl_ptr::Ptr{Cvoid}, x::Ptr{Float64},
    n::Ref{Int64},f::Ref{Float64},g::Ptr{Float64}, lind::Ptr{Vector{Int32}}, lbnd::Ptr{Float64},
    uind::Ptr{Vector{Int32}}, ubnd::Ptr{Float64},
    m::Ref{Int64}, p::Ref{Int64}, lambda::Ptr{Float64},
    jnnzmax::Ref{Float64}, hlnnzmax::Ref{Float64}, epsfeas::Ref{Float64},
    epscompl::Ref{Float64},epsopt::Ref{Float64}, rhoauto::Ref{Int32},
    rhoini::Ref{Float64},scale::Ref{Int32},extallowed::Ref{Int32},corrin::Ref{Int32},
    inform::Ref{Int64}
    )::Cvoid

  Libdl.dlclose(lib)
end