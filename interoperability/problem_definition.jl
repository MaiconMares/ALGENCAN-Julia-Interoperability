"Defines parameters and evaluation functions to optimization problems and export them to be used in ALGENCAN."
module ProblemDefinition
  export MyDataPtr, evalf!, evalc!, evalg!, evalj!, evalhl!, problem_params

  mutable struct MyDataPtr
    counters::NTuple{5, Int32}

    MyDataPtr(counters) = new(counters)
  end
  
  """
    problem_params()::Tuple

    Returns a Tuple of parameters for a CUTEst model with types expected by ALGENCAN.
  """
  function problem_params()::Tuple
    # Number of variables
    n::Int32 = 2

    # Initial guess and bound constraints

    x = Vector{Float64}([0, 0])

    # Stores the objective function
    f::Float64 = 0.0

    lind = Vector{Int32}(zeros(n))
    lbnd = Vector{Float64}(zeros(n))

    uind = Vector{Int32}(zeros(n))
    ubnd = Vector{Float64}(zeros(n))

    # Number of equality (m) and inequality (p) constraints
    m = 0
    p = 3

    # Initial guess for the Lagrange multipliers

    lambda = Vector{Float64}(zeros(m + p))

    # Number of entries in the Jacobian of the constraints

    jnnzmax::Int32 = 2 * (m+p)

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
    rhoini::Float64 = 6.9250264749392183e-310

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

    pdata::MyDataPtr = MyDataPtr((0, 0, 0, 0, 0))

    return x, n, f, lind, lbnd, uind, ubnd, m, p, lambda, jnnzmax, hlnnzmax, epsfeas, epscompl, epsopt,
    rhoauto, rhoini, scale, extallowed, corrin, inform, maxoutit, pdata
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
  function evalf!(n::Int32, x::Ptr{Float64}, f::Ptr{Float64}, inform::Int32, pdataptr::Ptr{MyDataPtr} = nothing)::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)

    d1::Float64 = sqrt((-5.0 - x_wrap[1])^2.0 + (10.0 - x_wrap[2])^2.0)
    d2::Float64 = sqrt((2.0 - x_wrap[1])^2.0 + (1.0 - x_wrap[2])^2.0)
    d3::Float64 = sqrt((10.0 - x_wrap[1])^2.0 + (5.0 - x_wrap[2])^2.0)

    temp::Float64 = d1 + d2 + d3

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
  function evalg!(n::Int32, x::Ptr{Float64}, g::Ptr{Float64}, inform::Int32, pdataptr::Ptr{MyDataPtr} = nothing)::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)
    g_wrap::Vector{Float64} = unsafe_wrap(Array, g, n)

    dfdx1::Float64 = (x_wrap[1] + 5.0) / sqrt(x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])
    dfdx2::Float64 = (x_wrap[1] - 2.0) / sqrt(x_wrap[1]^2.0 - 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])
    dfdx3::Float64 = (x_wrap[1] - 10.0) / sqrt(x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])

    dfdy1::Float64 = (x_wrap[2] - 10.0) / sqrt(x_wrap[2]^2.0 - 20.0 * x_wrap[2] + x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 125.0)
    dfdy2::Float64 = (x_wrap[2] - 1.0) / sqrt(x_wrap[2]^2.0 - 2.0 * x_wrap[2] + x_wrap[1]^2.0 + 5.0 - 4.0 * x_wrap[1])
    dfdy3::Float64 = (x_wrap[2] - 5.0) / sqrt(x_wrap[2]^2.0 - 10.0 * x_wrap[2] + x_wrap[1]^2.0 + 125.0 - 20.0 * x_wrap[1])


    g_wrap[1] = (dfdx1 + dfdx2 + dfdx3)
    g_wrap[2] = (dfdy1 + dfdy2 + dfdy3)

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
    n::Int32, x::Ptr{Float64}, m::Int32, p::Int32, c::Ptr{Float64}, inform::Int32, pdataptr::Ptr{MyDataPtr} = nothing,
  )::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr 
    if ((m + p) == 0)
      return nothing
    end

    x_wrap::Vector{Float64} = unsafe_wrap(Array, x, n)
    c_wrap::Vector{Float64} = unsafe_wrap(Array, c, (m + p))

    c_wrap[1] = sqrt((-5.0 - x_wrap[1])^2.0 + (10.0 - x_wrap[2])^2.0) - 10
    c_wrap[2] = sqrt((2.0 - x_wrap[1])^2.0 + (1.0 - x_wrap[2])^2.0) - 10
    c_wrap[3] = sqrt((10.0 - x_wrap[1])^2.0 + (5.0 - x_wrap[2])^2.0) - 10

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
  function evalj!(n::Int32, x::Ptr{Float64}, m::Int32, p::Int32, ind::Ptr{Int32},
    sorted::Ptr{Int32}, jsta::Ptr{Int32}, jlen::Ptr{Int32}, lim::Int32,
    jvar::Ptr{Int32}, jval::Ptr{Float64}, inform::Int32, pdataptr::Ptr{MyDataPtr} = nothing,
  )::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    if ((m + p) == 0)
      return nothing
    end

    x_wrap = unsafe_wrap(Array, x, n)
    ind_wrap = unsafe_wrap(Array, ind, m + p)
    sorted_wrap = unsafe_wrap(Array, sorted, m + p)
    jsta_wrap = unsafe_wrap(Array, jsta, m + p)
    jlen_wrap = unsafe_wrap(Array, jlen, m + p)
    jvar_wrap = unsafe_wrap(Array, jvar, lim)
    jval_wrap = unsafe_wrap(Array, jval, lim)

    error_code::Int32 = -94

    dg1dx::Float64 = (x_wrap[1] + 5.0) / sqrt(x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])
    dg1dy::Float64 = (x_wrap[2] - 10.0) / sqrt(x_wrap[2]^2.0 - 20.0 * x_wrap[2] + x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 125.0)

    dg2dx::Float64 = (x_wrap[1] - 2.0) / sqrt(x_wrap[1]^2.0 - 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])
    dg2dy::Float64 = (x_wrap[2] - 1.0) / sqrt(x_wrap[2]^2.0 - 2.0 * x_wrap[2] + x_wrap[1]^2.0 + 5.0 - 4.0 * x_wrap[1])

    dg3dx::Float64 = (x_wrap[1] - 10.0) / sqrt(x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])
    dg3dy::Float64 = (x_wrap[2] - 5.0) / sqrt(x_wrap[2]^2.0 - 10.0 * x_wrap[2] + x_wrap[1]^2.0 + 125.0 - 20.0 * x_wrap[1])

    sorted_wrap[1:m+p] .= 0

    if (ind_wrap[1] != 0)
      if (lim < n)
        unsafe_store!(inform, error_code)
        return
      end

      jval_wrap[1] = dg1dx
      jval_wrap[2] = dg1dy

      jsta_wrap[1] = 1
      jlen_wrap[1] = n

      jvar_wrap[1] = 1
      jvar_wrap[2] = 2

      sorted_wrap[1] = 1
    end

    if (ind_wrap[2] != 0)
      if (lim < n)
        unsafe_store!(inform, error_code)
        return
      end

      jval_wrap[3] = dg2dx
      jval_wrap[4] = dg2dy

      jsta_wrap[2] = 3
      jlen_wrap[2] = n

      jvar_wrap[3] = 1
      jvar_wrap[4] = 2

      sorted_wrap[2] = 1
    end

    if (ind_wrap[3] != 0)
      if (lim < n)
        unsafe_store!(inform, error_code)
        return
      end

      jval_wrap[5] = dg3dx
      jval_wrap[6] = dg3dy

      jsta_wrap[3] = 1
      jlen_wrap[3] = n

      jvar_wrap[5] = 1
      jvar_wrap[6] = 2

      sorted_wrap[3] = 1
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
    n::Int32, x::Ptr{Float64}, m::Int32, p::Int32, lambda::Ptr{Float64}, lim::Int32,
    inclf::Int32, hlnnz::Ptr{Int32}, hlrow::Ptr{Int32}, hlcol::Ptr{Int32}, hlval::Ptr{Float64},
    inform::Ptr{Int32}, pdataptr::Ptr{MyDataPtr} = nothing,
  )::Nothing
    # I should define an equivalent call to c_f_pointer(pdataptr,pdata) and an equivalent structure to pdata and pdataptr
    x_wrap = unsafe_wrap(Array, x, n)
    hlcol_wrap = unsafe_wrap(Array, hlcol, lim)
    hlrow_wrap = unsafe_wrap(Array, hlrow, lim)
    hlval_wrap = unsafe_wrap(Array, hlval, lim)
    lambda_wrap = unsafe_wrap(Array, lambda, m + p)
    inform_wrap = unsafe_wrap(Array, inform, 1)
    error_code::Int32 = -95

    hlnnz_temp = 0

    println("x_wrap = $x_wrap")

    # If .not. inclf then the Hessian of the objective function must not be included

    if (Bool(inclf))
      if (hlnnz_temp + 2 > lim)
        unsafe_store!(inform, error_code)
        return
      end

      hlnnz_temp = hlnnz_temp + 1

      hlrow_wrap[hlnnz_temp] = 1
      hlcol_wrap[hlnnz_temp] = 1

      d2fdxx_num1::Float64 = (x_wrap[2] - 20.0 * x_wrap[2] + 100.0)
      d2fdxx_num2::Float64 = (x_wrap[2] - 2.0 * x_wrap[2] + 1.0)
      d2fdxx_num3::Float64 = (x_wrap[2] - 10.0 * x_wrap[2] + 25.0)

      d2fdxx_term1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])
      d2fdxx_term2::Float64 = (x_wrap[2]^2.0 - 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])
      d2fdxx_term3::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])
      hlval_wrap[hlnnz_temp] = d2fdxx_num1 / d2fdxx_term1 * sqrt(d2fdxx_term1) + d2fdxx_num2 / d2fdxx_term2 * sqrt(d2fdxx_term2) + d2fdxx_num3 / d2fdxx_term3 * sqrt(d2fdxx_term3)

      hlnnz_temp = hlnnz_temp + 1

      hlrow_wrap[hlnnz_temp] = 1
      hlcol_wrap[hlnnz_temp] = 1

      d2fdxy_num1::Float64 = (x_wrap[2] - 10.0) * (x_wrap[1] + 5.0)
      d2fdxy_num2::Float64 = (x_wrap[2] - 1.0) * (x_wrap[1] - 2.0)
      d2fdxy_num3::Float64 = (x_wrap[2] - 5.0) * (x_wrap[1] - 10.0)

      d2fdxy_den1::Float64 = (x_wrap[2]^2.0 - 20.0 * x_wrap[2] + x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 125.0)^3
      d2fdxy_den2::Float64 = (x_wrap[2] - 20.0 * x_wrap[2] + x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 125.0)^3
      d2fdxy_den3::Float64 = (x_wrap[2]^2.0 - 10.0 * x_wrap[2] + x_wrap[1]^2.0 + 125.0 - 20.0 * x_wrap[1])^3

      hlval_wrap[hlnnz_temp] = -1.0 * d2fdxy_num1 / sqrt(d2fdxy_den1) - 1.0 * d2fdxy_num2 / sqrt(d2fdxy_den2) - 1.0 * d2fdxy_num3 / sqrt(d2fdxy_den3)

      hlnnz_temp = hlnnz_temp + 1

      hlrow_wrap[hlnnz_temp] = 2
      hlcol_wrap[hlnnz_temp] = 2

      d2fdyy_num1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 25.0)
      d2fdyy_num2::Float64 = (x_wrap[1]^2.0 - 4.0 * x_wrap[1] + 4.0)
      d2fdyy_num3::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + 100.0)

      d2fdyy_term1::Float64 = (x_wrap[2]^2.0 - 20.0 * x_wrap[2] + x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 125.0)
      d2fdyy_term2::Float64 = (x_wrap[2]^2.0 - 2.0 * x_wrap[2] + x_wrap[1]^2.0 + 5.0 - 4.0 * x_wrap[1])
      d2fdyy_term3::Float64 = (x_wrap[2]^2.0 - 10.0 * x_wrap[2] + x_wrap[1]^2.0 + 125.0 - 20.0 * x_wrap[1])

      hlval_wrap[hlnnz_temp] = d2fdyy_num1 / d2fdyy_term1 * sqrt(d2fdyy_term1) + d2fdyy_num2 / d2fdyy_term2 * sqrt(d2fdyy_term2) + d2fdyy_num3 / d2fdyy_term3 * sqrt(d2fdyy_term3)

    end

    # Note that entries of the Hessian of the Lagrangian can be
    # repeated. If this is case, them sum of repeated entrances is
    # considered. This feature simplifies the construction of the
    # Hessian of the Lagrangian.

    if (hlnnz_temp + 1 > lim)
      unsafe_store!(inform, error_code)
      return
    end

    # second order gradient of g1 restriction

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 1
    hlcol_wrap[hlnnz_temp] = 1

    d2g1dxx_num1::Float64 = (x_wrap[2] - 20.0 * x_wrap[2] + 100.0)

    d2g1dxx_term1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[1] * (d2g1dxx_num1 / d2g1dxx_term1 * sqrt(d2g1dxx_term1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 1

    d2g1dxy_num1::Float64 = (x_wrap[2] - 10.0) * (x_wrap[1] + 5.0)

    d2g1dxy_den1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])^3

    hlval_wrap[hlnnz_temp] = lambda_wrap[1] * (-1.0 * d2g1dxy_num1 / sqrt(d2g1dxy_den1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 2

    d2g1dyy_num1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + 25.0)

    d2g1dyy_term1::Float64 = (x_wrap[1]^2.0 + 10.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 20.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[1] * (d2g1dyy_num1 / d2g1dyy_term1 * sqrt(d2g1dyy_term1))

    # second order gradient of g2 restriction

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 1
    hlcol_wrap[hlnnz_temp] = 1

    d2g2dxx_num1::Float64 = (x_wrap[2] - 2.0 * x_wrap[2] + 1.0)

    d2g2dxx_term1::Float64 = (x_wrap[1]^2.0 - 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[2] * (d2g2dxx_num1 / d2g2dxx_term1 * sqrt(d2g2dxx_term1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 1

    d2g2dxy_num1::Float64 = (x_wrap[2] - 1.0) * (x_wrap[1] - 2.0)

    d2g2dxy_den1::Float64 = (x_wrap[1]^2 - 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])^3

    hlval_wrap[hlnnz_temp] = lambda_wrap[2] * (-1.0 * d2g2dxy_num1 / sqrt(d2g2dxy_den1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 2

    d2g2dyy_num1::Float64 = (x_wrap[1]^2.0 - 4.0 * x_wrap[1] + 4.0)

    d2g2dyy_term1::Float64 = (x_wrap[1]^2.0 + 4.0 * x_wrap[1] + x_wrap[2]^2.0 + 5.0 - 2.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[2] * (d2g2dyy_num1 / d2g2dyy_term1 * sqrt(d2g2dyy_term1))

    # second order gradient of g3 restriction

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 1
    hlcol_wrap[hlnnz_temp] = 1

    d2g3dxx_num1::Float64 = (x_wrap[2] - 10.0 * x_wrap[2] + 25.0)

    d2g3dxx_term1::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[3] * (d2g3dxx_num1 / d2g3dxx_term1 * sqrt(d2g3dxx_term1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 1

    d2g3dxy_num1::Float64 = (x_wrap[2] - 5.0) * (x_wrap[1] - 10.0)

    d2g3dxy_den1::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])^3

    hlval_wrap[hlnnz_temp] = lambda_wrap[3] * (-1.0 * d2g3dxy_num1 / sqrt(d2g3dxy_den1))

    hlnnz_temp = hlnnz_temp + 1

    hlrow_wrap[hlnnz_temp] = 2
    hlcol_wrap[hlnnz_temp] = 2

    d2g3dyy_num1::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + 100.0)

    d2g3dyy_term1::Float64 = (x_wrap[1]^2.0 - 20.0 * x_wrap[1] + x_wrap[2]^2.0 + 125.0 - 10.0 * x_wrap[2])
    hlval_wrap[hlnnz_temp] = lambda_wrap[3] * (d2g3dyy_num1 / d2g3dyy_term1 * sqrt(d2g3dyy_term1))

    unsafe_store!(hlnnz, hlnnz_temp)

    nothing
  end
end
