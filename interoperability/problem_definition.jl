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

    # Stores the objective function
    f::Float64 = 0.0

    g = Vector{Float64}(zeros(n))

    lind = Vector{Int32}(ones(n))
    lbnd = Vector{Float64}(zeros(n))

    uind = Vector{Int32}(zeros(n))
    ubnd = Vector{Float64}(zeros(n))

    # Number of equality (m) and inequality (p) constraints
    m::Int64 = 1
    p::Int64 = 1

    # Defines which constraints should be included in jacobian matrix
    ind = Vector{Int32}(ones(m+p))

    # Initial guess for the Lagrange multipliers

    lambda = Vector{Float64}(zeros(m+p))

    # Number of entries in the Jacobian of the constraints

    jnnzmax::Int64 = 2 * n

    # This should be the number of entries in the Hessian of the
    # Lagrangian. But, in fact, some extra space is need (to store the
    # Hessian of the Augmented Lagrangian, whose size is hard to
    # predict, and/or to store the Jacobian of the KKT system). Thus,
    # declare it as large as possible.

    hlnnzmax::Int64 = typemax(Int64)

    # Feasibility, complementarity, and optimality tolerances

    epsfeas::Float64  = 1.0e-08
    epscompl::Float64 = 1.0e-08
    epsopt::Float64   = 1.0e-08

    maxoutit::Int64 = 50

    # rhoauto means that Algencan will automatically set the initial
    # value of the penalty parameter. If you set rhoauto = .false. then
    # you must set rhoini below with a meaningful value.

    rhoauto::Int32 = 1
    rhoini::Float64 = 0.0

    if !Bool(rhoauto)
      rhoini = 1.0e-08
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

    inform::Int64 = 0

    return x,n,f,g,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin,inform,ind
  end

  function evalf!(n::Int64,x,f::Ptr{Float64},inform::Ptr{Int64})::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
      x_wrap = unsafe_wrap(Array, x, n)

      temp = ( x_wrap[1] + 3.0 * x_wrap[2] + x_wrap[3] )^(2.0) + 4.0 * (x_wrap[1] - x_wrap[2])^(2.0)
      unsafe_store!(f, temp)

      nothing
  end

  function evalg!(n::Int64,x,g,inform::Int64)::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
      x_wrap = unsafe_wrap(Array, x, n)
      g_wrap = unsafe_wrap(Array, g, n)

      t1::Float64 = x_wrap[1] + 3.0 * x_wrap[2] + x_wrap[3]
      t2::Float64 = x_wrap[1] - x_wrap[2]
      g_wrap[1] = 2.0 * t1 + 8.0 * t2
      g_wrap[2] = 6.0 * t1 - 8.0 * t2
      g_wrap[3] = 2.0 * t1

      nothing
  end

  function evalc!(
    n::Int64,x,m::Int64,p::Int64,c,inform::Ptr{Int64}
    )::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
      x_wrap = unsafe_wrap(Array, x, n)
      c_wrap = unsafe_wrap(Array, c, m+p)

      c_wrap[1] = 1.0 - x_wrap[1] - x_wrap[2] - x_wrap[3]
      c_wrap[2] = - 6.0 * x_wrap[2] - 4.0 * x_wrap[3] + (x_wrap[1]^3.0) + 3.0

      nothing
  end

  function evalj!(n::Int64,x,m::Int64,p::Int64,ind,
    sorted,jsta,jlen,lim::Int64,
    jvar,jval,inform::Int64,pdataptr::MyDataPtr=nothing
    )::Nothing
      # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
      x_wrap = unsafe_wrap(Array, x, n)
      ind_wrap = unsafe_wrap(Array, ind, m+p)
      sorted_wrap = unsafe_wrap(Array, sorted, m+p)
      jsta_wrap = unsafe_wrap(Array, jsta, m+p)
      jlen_wrap = unsafe_wrap(Array, jlen, m+p)
      jvar_wrap = unsafe_wrap(Array, jvar, lim)
      jval_wrap = unsafe_wrap(Array, jval, lim)

      nnz1 = 0
      nnz2 = 0
      if (Bool(ind_wrap[1]))
        if ( lim < n )
            inform = -94
            return
        end

        jsta_wrap[1] = 1
        jlen_wrap[1] = n
        nnz1 = nnz1 + 1

        for i in 1:n
          jvar_wrap[i] = i
        end

        jval_wrap[1:n] .= -1.0
        nnz2 = nnz2 + n

        sorted_wrap[1] = 1
      end

      if (Bool(ind_wrap[2]))
        if ( lim < n )
          inform = -94
          return
        end

        jsta_wrap[nnz1+1] = nnz2 + 1
        jlen_wrap[nnz1+1] = n

        for i in 1:n
          jvar_wrap[nnz2 + i] = i
        end
        
        jval_wrap[nnz2+1] = - 3.0 * (x_wrap[1]^2.0)
        jval_wrap[nnz2+2] = 6.0
        jval_wrap[nnz2+3] = 4.0

        sorted_wrap[2] = 1
      end

      nothing
  end

  function evalhl!(
    n::Int64,x,m::Int64,p::Int64,lambda,lim::Int64,
    inclf::Int32,hlnnz::Int64,hlrow,hlcol,hlval,
    inform::Int64,pdataptr::MyDataPtr=nothing
    )::Nothing
    hlnnzmax::Int64 = typemax(Int64)
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    x_wrap = unsafe_wrap(Array, x, n)
    hlrow_wrap = unsafe_wrap(Array, hlrow, hlnnzmax)
    hlcol_wrap = unsafe_wrap(Array, hlcol, hlnnzmax)
    hlval_wrap = unsafe_wrap(Array, hlval, hlnnzmax)
    lambda_wrap = unsafe_wrap(Array, lambda, m+p)

    hlnnz = 0

    # If .not. inclf then the Hessian of the objective function must not be included

    if (Bool(inclf))
      if ( hlnnz + 2 > lim )
        inform = -95
        return
      end
      
      hlrow_wrap[1]= 1
      hlrow_wrap[2]= 2
      hlrow_wrap[3]= 2
      hlrow_wrap[4]= 3
      hlrow_wrap[5]= 3
      hlrow_wrap[6]= 3

      hlcol_wrap[1]= 1
      hlcol_wrap[2]= 1
      hlcol_wrap[3]= 2
      hlcol_wrap[4]= 1
      hlcol_wrap[5]= 2
      hlcol_wrap[6]= 3

      hlval_wrap[1]= 10.0
      hlval_wrap[2]= -2.0
      hlval_wrap[3]= 26.0
      hlval_wrap[4]= 2.0
      hlval_wrap[5]= 6.0
      hlval_wrap[6]= 2.0

      hlnnz = 6

    end

    # Note that entries of the Hessian of the Lagrangian can be
    # repeated. If this is case, them sum of repeated entrances is
    # considered. This feature simplifies the construction of the
    # Hessian of the Lagrangian.

    if ( hlnnz + 1 > lim )
      inform = -95
      return
    end

    hlrow_wrap[hlnnz+1] = 1
    hlcol_wrap[hlnnz+1] = 1
    hlval_wrap[hlnnz+1] = lambda_wrap[2] * ( - 6.0 * x_wrap[1] )

    nothing
  end
end