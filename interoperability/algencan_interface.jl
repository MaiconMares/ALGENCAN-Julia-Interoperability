#include("./problem_definition.jl")

#using .ProblemDefinition

using Libdl

curr_dir = Vector{String}(["./"])
lib_path = Libdl.find_library("libalgencanma.so", curr_dir)

function run_algencan()
  lib = Libdl.dlopen(lib_path, RTLD_NOW|RTLD_GLOBAL)
  evalf_ptr = @cfunction(evalf!, Nothing, (Ref{Int32},Ptr{Float64},Ptr{Float64},Ref{Int32},Ptr{MyDataPtr}))
  evalg_ptr = @cfunction(evalg!, Nothing, (Ref{Int32},Ptr{Float64},Ptr{Float64},Ref{Int32},Ptr{MyDataPtr}))
  evalc_ptr = @cfunction(evalc!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Float64},Ref{Int32}, Ptr{MyDataPtr}
  ))

  evalj_ptr = @cfunction(evalj!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Int32},
    Ptr{Int32},Ptr{Int32},Ptr{Int32},Ref{Int32},
    Ptr{Int32},Ptr{Float64},Ref{Int32},Ptr{MyDataPtr}
  ))
    
  evalhl_ptr = @cfunction(evalhl!, Nothing, (
    Ref{Int32},Ptr{Float64},Ref{Int32},Ref{Int32},Ptr{Float64},Ref{Int32},
    Ref{Int32},Ptr{Int32},Ptr{Int32},Ptr{Int32},Ptr{Float64},
    Ptr{Int32},Ptr{MyDataPtr}
  ))

  x,n,f,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin,inform,maxoutit,pdata = problem_params()

  @ccall lib_path.__algencanma_MOD_init(
    evalf_ptr::Ptr{Cvoid}, evalg_ptr::Ptr{Cvoid},
    evalc_ptr::Ptr{Cvoid}, evalj_ptr::Ptr{Cvoid},
    evalhl_ptr::Ptr{Cvoid}, x::Ptr{Float64},
    n::Ref{Int32},f::Ref{Float64},lind::Ptr{Int32}, lbnd::Ptr{Float64},
    uind::Ptr{Int32}, ubnd::Ptr{Float64},
    m::Ref{Int32}, p::Ref{Int32}, lambda::Ptr{Float64},
    jnnzmax::Ref{Int32}, hlnnzmax::Ref{Int32}, epsfeas::Ref{Float64},
    epscompl::Ref{Float64},epsopt::Ref{Float64}, rhoauto::Ref{Int32},
    rhoini::Ref{Float64},scale::Ref{Int32},extallowed::Ref{Int32},corrin::Ref{Int32},
    inform::Ref{Int32},maxoutit::Ref{Int32},pdata::Ref{MyDataPtr}
    )::Cvoid

  Libdl.dlclose(lib)
end