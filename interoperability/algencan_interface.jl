include("./problem_definition.jl")
using .ProblemDefinition
using Libdl

function run_algencan()
  lib = Libdl.dlopen("libalgencanma.so")
  evalf_ptr = @cfunction(evalf!, Nothing, (Ref{Int64},Ptr{Vector{Float64}},Ref{Float64},Ref{Int64},Ref{MyDataPtr}))
  evalg_ptr = @cfunction(evalg!, Nothing, (Ref{Int64},Ptr{Vector{Float64}},Ptr{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}))
  evalc_ptr = @cfunction(evalc!, Nothing, (
    Ref{Int64},Ptr{Vector{Float64}},Ref{Int64},Ref{Int64},Ptr{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}
  ))

  evalj_ptr = @cfunction(evalj!, Nothing, (
    Ref{Int64},Ptr{Vector{Float64}},Ref{Int64},Ref{Int64},Ref{Vector{Int32}},
    Ref{Vector{Int32}},Ref{Vector{Int64}},Ref{Vector{Int64}},Ref{Int64},
    Ref{Vector{Int64}},Ptr{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}
  ))

  evalhl_ptr = @cfunction(evalhl!, Nothing, (
    Ref{Int64},Ptr{Vector{Float64}},Ref{Int64},Ref{Int64},Ptr{Vector{Float64}},Ref{Int64},
    Ref{Int32},Ref{Int64},Ref{Vector{Int64}},Ref{Vector{Int64}},Ptr{Vector{Float64}},
    Ref{Int64},Ref{MyDataPtr}
  ))

  ccall(
    (:__algencanma_MOD_init, "libalgencanma.so"), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid},Ptr{Cvoid}, Ptr{Cvoid},Ptr{Cvoid}),
    evalf_ptr, evalg_ptr, evalc_ptr, evalj_ptr, evalhl_ptr)

  Libdl.dlclose(lib)
end