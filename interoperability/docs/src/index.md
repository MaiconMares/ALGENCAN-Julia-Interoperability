# This is MyProject Documentation

Welcome to the documentation page. 

!!! note "Julia interface with ALGENCAN"
    Documentation about Julia interface with ALGENCAN and Julia interface with CUTEst to test CUTEst problems over
    ALGENCAN.

```@docs
JuliaInterface4CUTEst
JuliaInterface4CUTEst.evalf!(n::Int32,x::Ptr{Float64},f::Ptr{Float64},inform::Int32,pdataptr::Ptr{JuliaInterface4CUTEst.MyDataPtr}=nothing)
JuliaInterface4CUTEst.evalc!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,c::Ptr{Float64},inform::Int32,pdataptr::Ptr{JuliaInterface4CUTEst.MyDataPtr}=nothing
)
JuliaInterface4CUTEst.evalg!(n::Int32,x::Ptr{Float64},g::Ptr{Float64},inform::Int32,pdataptr::Ptr{JuliaInterface4CUTEst.MyDataPtr}=nothing)
JuliaInterface4CUTEst.evalj!(n::Int32,x::Ptr{Float64},m::Int32,p::Int32,ind::Ptr{Int32},
    sorted::Ptr{Int32},jsta::Ptr{Int32},jlen::Ptr{Int32},lim::Int32,
    jvar::Ptr{Int32},jval::Ptr{Float64},inform::Int32,pdataptr::Ptr{JuliaInterface4CUTEst.MyDataPtr}=nothing
    )
JuliaInterface4CUTEst.evalhl!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,lambda::Ptr{Float64},lim::Int32,
    inclf::Int32,hlnnz::Ptr{Int32},hlrow::Ptr{Int32},hlcol::Ptr{Int32},hlval::Ptr{Float64},
    inform::Ptr{Int32},pdataptr::Ptr{JuliaInterface4CUTEst.MyDataPtr}=nothing
)
JuliaInterface4CUTEst.problem_params()
csc2csr(nlp_obj::CUTEstModel,x_wrap::Vector{Float64},val2::Vector{Float64},col::Vector{Int32},ind2::Vector{Int32},line_len::Vector{Int32})

JuliaInterface4CUTEstTest
run_tests()

ProblemDefinition
ProblemDefinition.evalf!(n::Int32,x::Ptr{Float64},f::Ptr{Float64},inform::Int32,pdataptr::Ptr{ProblemDefinition.MyDataPtr}=nothing)
ProblemDefinition.evalc!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,c::Ptr{Float64},inform::Int32,pdataptr::Ptr{ProblemDefinition.MyDataPtr}=nothing
)
ProblemDefinition.evalg!(n::Int32,x::Ptr{Float64},g::Ptr{Float64},inform::Int32,pdataptr::Ptr{ProblemDefinition.MyDataPtr}=nothing)
ProblemDefinition.evalj!(n::Int32,x::Ptr{Float64},m::Int32,p::Int32,ind::Ptr{Int32},
    sorted::Ptr{Int32},jsta::Ptr{Int32},jlen::Ptr{Int32},lim::Int32,
    jvar::Ptr{Int32},jval::Ptr{Float64},inform::Int32,pdataptr::Ptr{ProblemDefinition.MyDataPtr}=nothing
    )
ProblemDefinition.evalhl!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,lambda::Ptr{Float64},lim::Int32,
    inclf::Int32,hlnnz::Ptr{Int32},hlrow::Ptr{Int32},hlcol::Ptr{Int32},hlval::Ptr{Float64},
    inform::Ptr{Int32},pdataptr::Ptr{ProblemDefinition.MyDataPtr}=nothing
)
ProblemDefinition.problem_params()
```