mutable struct mydataptr
end

function evalf(n::Int64,x::Vector{Float64},f::Float64, inform::Int64, pdataptr::mydataptr)
end

function evalg(n::Int64,x::Vector{Float64}, g, inform::Int64, pdataptr::mydataptr)
end

function evalc(n::Int64,x::Vector{Float64}, m, p, c, inform::Int64, pdataptr::mydataptr)
end

function evalj(n::Int64,x::Vector{Float64},m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform::Int64,pdataptr::mydataptr)
end

function evalhl(n::Int64,x::Vector{Float64},m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform::Int64,pdataptr::mydataptr)
end