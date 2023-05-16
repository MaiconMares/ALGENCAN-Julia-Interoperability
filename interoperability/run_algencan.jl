include("./algencan_interface.jl")
using .AlgencanInterface
using Libdl

function run_algencan()
  lib = Libdl.dlopen("libalgencanma.so")
  evalf_ptr = @cfunction(evalf!, nothing, (Ref{Int64},Ref{Vector{Float64}},Ref{Float64},Ref{Int64},Ref{MyDataPtr}))
  evalg_ptr = @cfunction(evalg!, nothing, (Ref{Int64},Ref{Vector{Float64}},Ref{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}))
  evalc_ptr = @cfunction(evalc!, nothing, (
    Ref{Int64},Ref{Vector{Float64}},Ref{Int64},Ref{Int64},Ref{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}
  ))

  evalj_ptr = @cfunction(evalj!, nothing, (
    Ref{Int64},Ref{Vector{Float64}},Ref{Int64},Ref{Int64},Ref{Vector{Int32}},
    Ref{Vector{Int32}},Ref{Vector{Int64}},Ref{Vector{Int64}},Ref{Int64},
    Ref{Vector{Int64}},Ref{Vector{Float64}},Ref{Int64},Ref{MyDataPtr}
  ))

  evalhl_ptr = @cfunction(evalhl!, nothing, (
    Ref{Int64},Ref{Vector{Float64}},Ref{Int64},Ref{Int64},Ref{Vector{Float64}},Ref{Int64},
    Ref{Int32},Ref{Int64},Ref{Vector{Int64}},Ref{Vector{Int64}},Ref{Vector{Float64}},
    Ref{Int64},Ref{MyDataPtr}
  ))

  ccall(
    (:init, "libalgencanma.so"), Cvoid,
    (Ptr{Cvoid}, Ptr{Cvoid},Ptr{Cvoid}, Ptr{Cvoid},Ptr{Cvoid}),
    evalf_ptr, evalg_ptr, evalc_ptr, evalj_ptr, evalhl_ptr)

  Libdl.dlclose(lib)
end