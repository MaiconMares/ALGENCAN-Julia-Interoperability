module ProblemDefinition
  export MyDataPtr,evalf!,evalc!,evalg!,evalj!,evalhl!,problem_params

  mutable struct MyDataPtr
    counters::NTuple{5, Cint}

    MyDataPtr(counters) = new((0,0,0,0,0))
  end

  function problem_params()::Tuple
    # Number of variables
    n::Int64 = 3

    # Initial guess and bound constraints

    x = Vector{Float64}([0.1,0.7,0.2])

    lind = Vector{Int32}(ones(n))
    lbnd = Vector{Float64}(zeros(n))

    uind = Vector{Int32}(zeros(n))
    ubnd = Vector{Float64}(zeros(n))

    # Number of equality (m) and inequality (p) constraints
    m::Int64 = 1
    p::Int64 = 1

    # Initial guess for the Lagrange multipliers

    lambda = Vector{Float64}(zeros(m+p))

    # Number of entries in the Jacobian of the constraints

    jnnzmax::Int64 = 2 * n

    # This should be the number of entries in the Hessian of the
    # Lagrangian. But, in fact, some extra space is need (to store the
    # Hessian of the Augmented Lagrangian, whose size is hard to
    # predict, and/or to store the Jacobian of the KKT system). Thus,
    # declare it as large as possible.

    hlnnzmax = typemax(Int64)

    # Feasibility, complementarity, and optimality tolerances

    epsfeas  = 1.0e-08
    epscompl = 1.0e-08
    epsopt   = 1.0e-08

    maxoutit::Int64 = 50

    # rhoauto means that Algencan will automatically set the initial
    # value of the penalty parameter. If you set rhoauto = .false. then
    # you must set rhoini below with a meaningful value.

    rhoauto::Int32 = 1

    if !Bool(rhoauto)
      rhoini::Float64 = 1.0e-08
    end

    # scale = .true. means that you allow Algencan to automatically
    # scale the constraints. In any case, the feasibility tolerance
    # (epsfeas) will be always satisfied by the UNSCALED original
    # constraints.
    scale::Int32 = 0

    # extallowed = .true. means that you allow Gencan (the active-set
    # method used by Algencan to solve the bound-constrained
    # subproblems) to perform extrapolations. This strategy may use
    # extra evaluations of the objective function and the constraints
    # per iterations; but it uses to provide overal savings. You should
    # test both choices for the problem at hand.
    extallowed::Int32 = 1

    # corrin = .true. means that you allow the inertia of the
    # Jacobian of the KKT system to be corrected during the acceleration
    # process. You should test both choices for the problem at hand.
    corrin::Int32 = 0

    return x,n,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin
  end

  function evalf!(n::Int64,x::Vector{Float64},f::Float64,inform::Int64,pdataptr::MyDataPtr=nothing)::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
      println("evalf!() has been called!")

      f = ( x[1] + 3.0 * x[2] + x[3] )^(2.0) + 4.0 * (x[1] - x[2])^(2.0)

      println(f)
  end

  function evalg!(n::Int64,x::Vector{Float64},g::Vector{Float64},inform::Int64,pdataptr::MyDataPtr=nothing)::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      t1::Int64 = x[1] + 3.0 * x[2] + x[3]
      t2::Int64 = x[1] - x[2]
      g[1] = 2.0 * t1 + 8.0 * t2
      g[2] = 6.0 * t1 - 8.0 * t2
      g[3] = 2.0 * t1
  end

  function evalc!(
    n::Int64,x::Vector{Float64},m::Int64,p::Int64,c::Vector{Float64},inform::Int64,
    pdataptr::MyDataPtr=nothing
    )::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr

      c[1] = 1.0 - x[1] - x[2] - x[3]
      c[2] = - 6.0 * x[2] - 4.0 * x[3] + (x[1]^3.0) + 3.0
  end

  function evalj!(n::Int64,x::Vector{Float64},m::Int64,p::Int64,ind::Vector{Int32},
    sorted::Vector{Int32},jsta::Vector{Int64},jlen::Vector{Int64},lim::Int64,
    jvar::Vector{Int64},jval::Vector{Float64},inform::Int64,pdataptr::MyDataPtr=nothing
    )::Nothing
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

  function evalhl!(
    n::Int64,x::Vector{Float64},m::Int64,p::Int64,lambda::Vector{Float64},lim::Int64,
    inclf::Int32,hlnnz::Int64,hlrow::Vector{Int64},hlcol::Vector{Int64},hlval::Vector{Float64},
    inform::Int64,pdataptr::MyDataPtr=nothing
    )::Nothing
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