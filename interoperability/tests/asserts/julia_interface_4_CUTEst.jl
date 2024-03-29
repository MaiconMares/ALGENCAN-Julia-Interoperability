"A module intended to define parameters and functions for CUTEst models in a way that can be used with any model.
It exports parameters and evaluation functions to be used in ALGENCAN.
It uses NLPModels package to compute derivatives and another needed computations."
module JuliaInterface4CUTEst
  using NLPModels, CUTEst
  export MyDataPtr,evalf!,evalc!,evalg!,evalj!,evalhl!,problem_params,csc2csr, nlp

  mutable struct MyDataPtr
    counters::NTuple{5, Int32}

    MyDataPtr(counters) = new(counters)
  end

  nlp = ""

  """
    problem_params()::Tuple

    Returns a Tuple of parameters for a CUTEst model with types expected by ALGENCAN.
  """
  function problem_params()::Tuple
    # Number of variables
    n::Int32 = nlp.meta.nvar

    # Initial guess and bound constraints

    x = Vector{Float64}(nlp.meta.x0)

    # Stores the objective function
    f::Float64 = 0.0

    lind = Vector{Int32}(zeros(n))
    lbnd::Vector{Float64} = [-1.0e20 for _ in 1:n]

    uind = Vector{Int32}(zeros(n))
    ubnd::Vector{Float64} = [1.0e20 for _ in 1:n]

    lvar = nlp.meta.lvar
    uvar = nlp.meta.uvar
    for k in 1:length(lvar)
      if (lvar[k] != -Inf)
        lind[k] = 1
        lbnd[k] = lvar[k]
      end

      if (uvar[k] != Inf)
        uind[k] = 1
        ubnd[k] = uvar[k]
      end
    end

    # Number of equality (m) and inequality (p) constraints
    m = length(nlp.meta.jfix)
    p = length(nlp.meta.jlow) + length(nlp.meta.jupp) + length(nlp.meta.jrng)

    # Initial guess for the Lagrange multipliers

    lambda = Vector{Float64}(zeros(m+p))

    # Number of entries in the Jacobian of the constraints

    jnnzmax::Int32 = 2 * nlp.meta.nnzj

    # This should be the number of entries in the Hessian of the
    # Lagrangian. But, in fact, some extra space is need (to store the
    # Hessian of the Augmented Lagrangian, whose size is hard to
    # predict, and/or to store the Jacobian of the KKT system). Thus,
    # declare it as large as possible.

    hlnnzmax::Int32 = typemax(Int32) ÷ 10

    # Feasibility, complementarity, and optimality tolerances

    epsfeas::Float64  = 1.0e-08
    epscompl::Float64 = 1.0e-08
    epsopt::Float64   = 1.0e-08

    maxoutit::Int32 = 50

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

    inform::Int32 = 0

    pdata::MyDataPtr = MyDataPtr((0,0,0,0,0))

    return x,n,f,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,epsopt,
          rhoauto,rhoini,scale,extallowed,corrin,inform,maxoutit,pdata
  end

  """
    evalf!(n::Int32,x::Ptr{Float64},f::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing)::Nothing

    Computes objective function with parameters provided from ALGENCAN and stores it in the same 
    address of f::Ptr{Float64}.
    Parameters:
    n::Int32: number of variables in optimization problem (dimension of vector x).
    x::Ptr{Float64}: vector that stores variables of optimization problem.
    f::Ptr{Float64}: stores value of objective function.
    inform::Int32: holds error code if any error occurrs during solving.
    pdataptr::Ptr{MyDataPtr}: computes number of calls to a function.
  """
  function evalf!(n::Int32,x::Ptr{Float64},f::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing)::Nothing
    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)

    temp::Float64 = convert(Float64, obj(nlp, x_wrap))
    
    unsafe_store!(f, temp)

    nothing
  end

  """
    evalg!(n::Int32,x::Ptr{Float64},g::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing)::Nothing

    Computes gradient of objective function with parameters provided from ALGENCAN and stores it in the same 
    address of g::Ptr{Float64}.
    Parameters:
    n::Int32: number of variables in optimization problem (dimension of vector x).
    x::Ptr{Float64}: vector that stores variables of optimization problem.
    g::Ptr{Float64}: stores gradient of objective function.
    inform::Int32: holds error code if any error occurrs during solving.
    pdataptr::Ptr{MyDataPtr}: computes number of calls to a function.
  """
  function evalg!(n::Int32,x::Ptr{Float64},g::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing)::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)
    g_wrap::Vector{Float64} = unsafe_wrap(Array, g, n)
    
    df = (x_arg::Vector{Float64}) -> convert(Vector{Float64},grad(nlp, x_arg))

    temp::Vector{Float64} = df(x_wrap)
    
    for j in 1:n
      g_wrap[j] = temp[j]
    end

    nothing
  end

  """
    evalc!(
      n::Int32,x::Ptr{Float64},m::Int32,p::Int32,c::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing

    Computes vector of restrictions of size m+p with parameters provided from ALGENCAN. Stores result in the same 
    address of c::Ptr{Float64}.
    Parameters:
    n::Int32: number of variables in optimization problem (dimension of vector x).
    x::Ptr{Float64}: vector that stores variables of optimization problem.
    m::Int32: number of equality constraints.
    p::Int32: number of inequality constraints.
    c::Ptr{Float64}: stores computed constraints, equality constraints comes first.
    inform::Int32: holds error code if any error occurrs during solving.
    pdataptr::Ptr{MyDataPtr}: computes number of calls to a function.
  """
  function evalc!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,c::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing
    if ((m+p) == 0)
      return nothing
    end

    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)
    c_wrap::Vector{Float64} = unsafe_wrap(Array, c, (m+p))

    temp::Vector{Float64} = convert(Vector{Float64}, cons(nlp, x_wrap))

    low_constr_idx = nlp.meta.jlow
    for j in 1:(m+p)
      if j in low_constr_idx
        c_wrap[j] = temp[j]*-1.0
      else
        c_wrap[j] = temp[j]
      end
    end

    nothing
  end

  """
    evalj!(n::Int32,x::Ptr{Float64},m::Int32,p::Int32,ind::Ptr{Int32},
      sorted::Ptr{Int32},jsta::Ptr{Int32},jlen::Ptr{Int32},lim::Int32,
      jvar::Ptr{Int32},jval::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing

    Computes vectors that describes jacobian of restrictions in Compressed Sparse Row format. Computed vectors are
    stored in the same addresses of jsta::Ptr{Int32},jlen::Ptr{Int32},jvar::Ptr{Int32} and jval::Ptr{Float64}.
    Parameters:
    n::Int32: number of variables in optimization problem (dimension of vector x).
    x::Ptr{Float64}: vector that stores variables of optimization problem.
    m::Int32: number of equality constraints.
    p::Int32: number of inequality constraints.
    ind::Ptr{Int32}: boolean that indicates which constraints should be included in jacobian of constraints.
    sorted::Ptr{Int32}: indicates if gradients are in ascending order.
    jsta::Ptr{Int32}: stores indexes from jval of elements that are first elements of its respective line in original matrix.
    jlen::Ptr{Int32}: stores quantity of non zero elements in each line.
    lim::Int32: limit of elements in jvar and jval.
    jvar::Ptr{Int32}: stores columns indexes of elements from jval.
    jval::Ptr{Float64}: stores non zero elements.
    inform::Int32: holds error code if any error occurrs during solving.
    pdataptr::Ptr{MyDataPtr}: computes number of calls to a function.
  """
  function evalj!(n::Int32,x::Ptr{Float64},m::Int32,p::Int32,ind::Ptr{Int32},
    sorted::Ptr{Int32},jsta::Ptr{Int32},jlen::Ptr{Int32},lim::Int32,
    jvar::Ptr{Int32},jval::Ptr{Float64},inform::Int32,pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing
    if ((m+p) == 0)
      return nothing
    end

    x_wrap = unsafe_wrap(Array, x, n)
    ind_wrap = unsafe_wrap(Array, ind, m+p)
    sorted_wrap = unsafe_wrap(Array, sorted, m+p)
    jsta_wrap = unsafe_wrap(Array, jsta, m+p)
    jlen_wrap = unsafe_wrap(Array, jlen, m+p)
    jvar_wrap = unsafe_wrap(Array, jvar, lim)
    jval_wrap = unsafe_wrap(Array, jval, lim)
    
    csc2csr(nlp,x_wrap,jval_wrap,jvar_wrap,jsta_wrap,jlen_wrap)
    
    for i in 1:(m+p)
      sorted_wrap[i] = 0
      if (ind_wrap[i] != 0)
        if ( lim < n )
          error_code::Int32 = -94
          unsafe_store!(inform, error_code)
          return
        end
        
      end
    end
    
    low_constr_idx = nlp.meta.jlow
    if (length(low_constr_idx) > 0)
      for i in (m*n)+1:length(jval_wrap)
        jval_wrap[i] = jval_wrap[i] * -1.0
      end
    end

    nothing
  end

  """
    evalhl!(
      n::Int32,x::Ptr{Float64},m::Int32,p::Int32,lambda::Ptr{Float64},lim::Int32,
      inclf::Int32,hlnnz::Ptr{Int32},hlrow::Ptr{Int32},hlcol::Ptr{Int32},hlval::Ptr{Float64},
      inform::Ptr{Int32},pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing

    Computes vectors that describes hessian of Lagrangian function in Coordinate format. Computed vectors are
    stored in the same addresses of hlrow::Ptr{Int32},hlcol::Ptr{Int32} and hlval::Ptr{Float64}.
    hlrow stores indexes of line elements and hlcol their columns indexes. hlval stores non zero elements from
    original matrix.
    Parameters:
    n::Int32: number of variables in optimization problem (dimension of vector x)
    x::Ptr{Float64}: vector that stores variables of optimization problem.
    m::Int32: number of equality constraints.
    p::Int32: number of inequality constraints.
    lambda::Ptr{Float64}: stores lagrangian multipliers
    lim::Int32: limit of elements in hlrow, hlcol and hlval.
    inclf::Int32: indicates if objective function should be included in hessian.
    hlnnz::Ptr{Int32}: number of non-zero elements in hessian.
    hlrow::Ptr{Int32}: stores line indexes of elements in hlval.
    hlcol::Ptr{Int32}: stores column indexes of elements in hlval.
    hlval::Ptr{Float64}: stores non-zero elements of hessian sparse matrix.
    inform::Int32: holds error code if any error occurrs during solving.
    pdataptr::Ptr{MyDataPtr}: computes number of calls to a function.
  """
  function evalhl!(
    n::Int32,x::Ptr{Float64},m::Int32,p::Int32,lambda::Ptr{Float64},lim::Int32,
    inclf::Int32,hlnnz::Ptr{Int32},hlrow::Ptr{Int32},hlcol::Ptr{Int32},hlval::Ptr{Float64},
    inform::Ptr{Int32},pdataptr::Ptr{MyDataPtr}=nothing
    )::Nothing
    x_wrap = unsafe_wrap(Array, x, n)
    hlcol_wrap = unsafe_wrap(Array, hlcol, lim)
    hlval_wrap = unsafe_wrap(Array, hlval, lim)
    lambda_wrap = unsafe_wrap(Array, lambda, m+p)
    inform_wrap = unsafe_wrap(Array, inform, 1)
    error_code::Int32 = -95

    hlnnz_temp = 0

    # If .not. inclf then the Hessian of the objective function must not be included
    
    if (Bool(inclf))
      if ( hlnnz_temp + 2 > lim )
        unsafe_store!(inform, error_code)
        return
      end

      low_constr_idx = nlp.meta.jlow
      if (length(low_constr_idx) > 0)
        for i in low_constr_idx
          lambda_wrap[i] = lambda_wrap[i] * -1.0
        end
      end

      hlval_temp = hess_coord(nlp, x_wrap, lambda_wrap)
      hlrow_temp, hlcol_temp = hess_structure(nlp)

      for val in 1:length(hlval_temp)
        unsafe_store!(hlval, hlval_temp[val], val)
      end

      for row in 1:length(hlrow_temp)
        unsafe_store!(hlcol, hlrow_temp[row], row)
      end

      for col in 1:length(hlcol_temp)
        unsafe_store!(hlrow, hlcol_temp[col], col)
      end

      hlnnz_temp = nlp.meta.nnzh
      unsafe_store!(hlnnz, hlnnz_temp)
    end

    # Note that entries of the Hessian of the Lagrangian can be
    # repeated. If this is case, them sum of repeated entrances is
    # considered. This feature simplifies the construction of the
    # Hessian of the Lagrangian.

    if (hlnnz_temp + 1 > lim)
      unsafe_store!(inform, error_code)
      return
    end

    nothing
  end

  """
    csc2csr(
      nlp_obj::CUTEstModel,x_wrap::Vector{Float64},val2::Vector{Float64},col::Vector{Int32},
      ind2::Vector{Int32},line_len::Vector{Int32}
    )::Nothing

    Computes jacobian of constraints and convert it to Compressed Sparse Row format.
    Parameters:
    nlp_obj::CUTEstModel: stores a CUTEst object that describes optimization problem.
    x_wrap::Vector{Float64}: vector that stores variables of optimization problem
    val2::Vector{Float64}: stores non zero elements in CSC format.
    col::Vector{Int32}: stores columns indexes of elements from jval.
    ind2::Vector{Int32}: stores indexes from jval of elements that are first elements 
    of its respective line in original matrix.
    line_len::Vector{Int32}: stores quantity of non zero elements in each line.
  """
  function csc2csr(
    nlp_obj::CUTEstModel,x_wrap::Vector{Float64},val2::Vector{Float64},col::Vector{Int32},
    ind2::Vector{Int32},line_len::Vector{Int32}
  )::Nothing
    jacobian = jac(nlp_obj, x_wrap)
    m::Int64 = jacobian.m
    n::Int64 = jacobian.n
    ind1::Vector{Int64} = jacobian.colptr
    val1::Vector{Float64} = jacobian.nzval
    row::Vector{Int64} = jacobian.rowval
    nnz::Int64 = length(jacobian.nzval)
    cnt = Vector{Int64}(zeros(m))

    # Cross row to check how many non zero elements there are in each line
    for i in 1:nnz
      cnt[row[i]] = cnt[row[i]] + 1
    end

    # Fill ind2 with indexes indicating where starts each line according to val1
    ind2[1] = 1
    for i in 2:m
      ind2[i] = ind2[i-1] + cnt[i-1]

      # Copy num of elements per line from cnt to jlen
      line_len[i-1] = cnt[i-1]

      # Copy ind2 to cnt
      cnt[i-1] = ind2[i-1]
    end
    line_len[m] = cnt[m]
    cnt[m] = ind2[m]

    # Fill col and val2 using cnt and row as guides
    for j in 1:n-1
      for i in ind1[j]:ind1[j+1]-1
        col[cnt[row[i]]] = j
        val2[cnt[row[i]]] = val1[i]
        cnt[row[i]] = cnt[row[i]] + 1
      end
    end

    # Cross the last column
    for j in ind1[n]:nnz
      col[cnt[row[j]]] = n
      val2[cnt[row[j]]] = val1[j]
      cnt[row[j]] = cnt[row[j]] + 1
    end

    nothing
  end
end