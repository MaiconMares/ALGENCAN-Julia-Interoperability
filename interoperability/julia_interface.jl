mutable struct mydataptr
end

function evalf(n::Int64,x::Vector{Float64},f::Float64,inform::Int64,pdataptr::mydataptr)
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

    f = ( x[1] + 3.0 * x[2] + x[3] )^(2.0) + 4.0 * (x[1] - x[2])^(2.0)
end

function evalg(n::Int64,x::Vector{Float64},g::Vector{Float64},inform::Int64,pdataptr::mydataptr)
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

    t1::Int64 = x[1] + 3.0 * x[2] + x[3]
    t2::Int64 = x[1] - x[2]
    g[1] = 2.0 * t1 + 8.0 * t2
    g[2] = 6.0 * t1 - 8.0 * t2
    g[3] = 2.0 * t1
end

function evalc(n::Int64,x::Vector{Float64},m::Int64,p::Int64,c::Vector{Float64},inform::Int64,pdataptr::mydataptr)
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

    c[1] = 1.0 - x[1] - x[2] - x[3]
    c[2] = - 6.0 * x[2] - 4.0 * x[3] + (x[1]^3.0) + 3.0
end

function evalj(n::Int64,x::Vector{Float64},m::Int64,p::Int64,ind::Vector{Int32},sorted::Vector{Int32},jsta::Vector{Int64},jlen::Vector{Int64},lim::Int64,jvar::Vector{Int64},jval::Vector{Float64},inform::Int64,pdataptr::mydataptr)
    nnz1::Int64 = 0
    nnz2::Int64 = 0

    if ( ind[1] )
       if ( lim < n )
          inform = -94
          return
       end

       jsta[1] = 1
       jlen[1] = n
       nnz1 = nnz1 + 1

       jvar = [i for i in 1:n]
       jval[1:n] = -1.0
       nnz2 = nnz2 + n

       sorted[1] = 1
    end

    if ( ind[2] )
       if ( lim .lt. n )
          inform = -94
          return
       end

       jsta(nnz1+1) = nnz2 + 1
       jlen(nnz1+1) = n

       jvar(nnz2+1:nnz2+n) = (/ (i,i=1,n) /)

       jval(nnz2+1) = - 3.0d0 * x[1] ** 2.0d0
       jval(nnz2+2) = 6.0d0
       jval(nnz2+3) = 4.0d0

       sorted[2] = 1
    end
end

function evalhl(n::Int64,x::Vector{Float64},m::Int64,p::Int64,lambda::Vector{Float64},lim::Int64,inclf::Int32,hlnnz::Int64,hlrow::Vector{Int64},hlcol::Vector{Int64},hlval::Vector{Float64},inform::Int64,pdataptr::mydataptr)
end