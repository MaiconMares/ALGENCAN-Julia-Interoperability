module AlgencanInterface
  export MyDataPtr, evalf, evalc, evalg, evalj, evalhl
  
  mutable struct MyDataPtr
    counters::NTuple{5, Cint}

    MyDataPtr(counters) = new((0,0,0,0,0))
  end

  function evalf(n::Int64,x::Vector{Float64},f::Float64,inform::Int64,pdataptr::MyDataPtr=nothing)
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      f = ( x[1] + 3.0 * x[2] + x[3] )^(2.0) + 4.0 * (x[1] - x[2])^(2.0)
  end

  function evalg(n::Int64,x::Vector{Float64},g::Vector{Float64},inform::Int64,pdataptr::MyDataPtr=nothing)
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      t1::Int64 = x[1] + 3.0 * x[2] + x[3]
      t2::Int64 = x[1] - x[2]
      g[1] = 2.0 * t1 + 8.0 * t2
      g[2] = 6.0 * t1 - 8.0 * t2
      g[3] = 2.0 * t1
  end

  function evalc(
    n::Int64,x::Vector{Float64},m::Int64,p::Int64,c::Vector{Float64},inform::Int64,
    pdataptr::MyDataPtr=nothing
    )
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      c[1] = 1.0 - x[1] - x[2] - x[3]
      c[2] = - 6.0 * x[2] - 4.0 * x[3] + (x[1]^3.0) + 3.0
  end

  function evalj(n::Int64,x::Vector{Float64},m::Int64,p::Int64,ind::Vector{Int32},
    sorted::Vector{Int32},jsta::Vector{Int64},jlen::Vector{Int64},lim::Int64,
    jvar::Vector{Int64},jval::Vector{Float64},inform::Int64,pdataptr::MyDataPtr=nothing
    )
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      if ( ind[1] )
        if ( lim < n )
            inform = -94
            return
        end

        jsta[1] = 1
        jlen[1] = n
        nnz1 = nnz1 + 1

        jvar = [i for i in 1:n]
        jval[1:n] .= -1.0

        sorted[1] = 1
      end

      if ( ind[2] )
        if ( lim < n )
            inform = -94
            return
        end

        push!(jsta, length(jval) + 1)
        push!(jlen, n)

        for i in 1:n
          push!(jvar, i)
        end

        push!(jval, - 3.0 * (x[1]^2.0))
        push!(jval, 6.0)
        push!(jval, 4.0)

        sorted[2] = 1
      end
  end

  function evalhl(
    n::Int64,x::Vector{Float64},m::Int64,p::Int64,lambda::Vector{Float64},lim::Int64,
    inclf::Int32,hlnnz::Int64,hlrow::Vector{Int64},hlcol::Vector{Int64},hlval::Vector{Float64},
    inform::Int64,pdataptr::MyDataPtr=nothing
    )
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

    hlnnz = 0

    # If .not. inclf then the Hessian of the objective function must not be included

    if ( inclf )
      if ( hlnnz + 2 > lim )
        inform = -95
        return
      end

      hlnnz = 6

      hlrow = [      1,      2,      2,     3,     3,     3 ]
      hlcol = [      1,      1,      2,     1,     2,     3 ]
      hlval = [ 10.0, -2.0, 26.0, 2.0, 6.0, 2.0 ]

    end

    # Note that entries of the Hessian of the Lagrangian can be
    # repeated. If this is case, them sum of repeated entrances is
    # considered. This feature simplifies the construction of the
    # Hessian of the Lagrangian.

    if ( hlnnz + 1 > lim )
      inform = -95
      return
    end

    hlnnz = hlnnz + 1

    push!(hlrow, 1)
    push!(hlcol, 1)
    push!(hlval, lambda[2] * ( - 6.0 * x[1] ))
  end
end