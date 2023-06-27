! *****************************************************************
! *****************************************************************

! This is a main program that calls Algencan to solve a simple
! problem. It is intended to be used as an (incomplete) example of
! usage. Algencan applies to problems of the form
!
! Minimize f(x)
!
! subject to
!
!   heq(x) = 0    heq : R^n \to R^m represents the equality constraints
!   hin(x) <= 0   hin : R^n \to R^p represents the inequality constraints
!   l <= x <= u   l and u \in R^n are the bound constraints.

! *****************************************************************
! *****************************************************************

module algencanma

  use bmalgencan, only: algencan
  use iso_c_binding, only: c_ptr, c_loc,c_f_pointer

  implicit none

  ! Re-define this type (pdata_type) anyway you want. Algencan
    ! receives a pointer to a 'structure of this type'. Algencan has no
    ! access to the structure. It simple passes the pointer back to the
    ! user defined subroutines evalf, evalg, evalc, evalj, and
    ! evalhl. So, this is a trade safe way (no common blocks) of passing
    ! to the user-provided routines any information related to the
    ! problem. In this example, it is only being used for the user to
    ! count by itself the number of calls to each routine.
  type :: pdata_type
    integer :: counters(5) = 0
  end type pdata_type

  abstract interface
    subroutine evalf(n,x,f,inform,pdataptr)
      import c_ptr

      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(inout) :: inform
      real(kind=8), intent(out) :: f
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)

      ! This routine must compute the objective function.
    end subroutine evalf

    subroutine evalg(n,x,g,inform,pdataptr)
      import c_ptr

      ! SCALAR ARGUMENTS
      integer, intent(in) :: n
      integer, intent(inout) :: inform
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: g(n)

      ! This routine must compute the gradient of the objective
      ! function.
    end subroutine evalg

    subroutine evalc(n,x,m,p,c,inform,pdataptr)
      import c_ptr

      ! SCALAR ARGUMENTS
      integer, intent(in) :: m,n,p
      integer, intent(inout) :: inform
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: x(n)
      real(kind=8), intent(out) :: c(m+p)

      ! This routine must compute all the m+p constraints.
    end subroutine evalc

    subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)
      import c_ptr

      ! SCALAR ARGUMENTS
      integer, intent(in) :: lim,m,n,p
      integer, intent(inout) :: inform
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      logical, intent(in) :: ind(m+p)
      real(kind=8), intent(in) :: x(n)
      logical, intent(out) :: sorted(m+p)
      integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
      real(kind=8), intent(out) :: jval(lim)

    end subroutine evalj

    subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,hlnnzmax,pdataptr)
      import c_ptr

      logical, intent(in) :: inclf
      integer, intent(in) :: m,n,lim,p,hlnnzmax
      integer, intent(out) :: hlnnz
      integer, intent(inout) :: inform
      type(c_ptr), optional, intent(in) :: pdataptr

      ! ARRAY ARGUMENTS
      real(kind=8), intent(in) :: lambda(m+p),x(n)
      integer, intent(out) :: hlrow(lim),hlcol(lim)
      real(kind=8), intent(out) :: hlval(lim)

    end subroutine evalhl
  end interface

contains

  subroutine init(user_evalf,user_evalg,user_evalc,user_evalj,user_evalhl, &
    x,n,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,  &
    epsopt,rhoauto,rhoini,scale,extallowed,corrin)

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(evalf) :: user_evalf
    procedure(evalg) :: user_evalg
    procedure(evalc) :: user_evalc
    procedure(evalj) :: user_evalj
    procedure(evalhl) :: user_evalhl

    ! PARAMETERS
    integer, intent(in) :: n, m, p
    logical, intent(in) :: corrin,extallowed,rhoauto,scale, lind(n),uind(n)
    real(kind=8), intent(inout) :: lbnd(n),ubnd(n),lambda(m+p),x(n), epsfeas,epscompl,epsopt,rhoini
    integer, intent(inout) :: hlnnzmax, jnnzmax

    ! LOCAL SCALARS

    integer :: allocerr,ierr,inform,istop,maxoutit,nwcalls,nwtotit,outiter,totiter
    real(kind=8) :: bdsvio,csupn,finish,nlpsupn,ssupn,start,f
    type(pdata_type), target :: pdata

    real(kind=8), allocatable :: c(:)

    if ( allocerr .ne. 0 ) then
      print *, 'Allocation error.'
      stop
    end if

    allocate(c(m+p),stat=allocerr)

    if ( allocerr .ne. 0 ) then
      print *, 'Allocation error.'
      stop
    end if

    call cpu_time(start)

    !call algencan(user_evalf,user_evalg,user_evalc,user_evalj,user_evalhl,jnnzmax,hlnnzmax, &
    !    n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
    !    scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
    !    outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

    call cpu_time(finish)

    call user_evalf(n,x,f,inform)

    print *, ''
    print *, 'Number of variables                                   = ',n
    print *, 'Number of equality constraints                        = ',m
    print *, 'Number of inequality constraints                      = ',p

    print *, ''
    print *, '(REPORTED BY SOLVER) istop                            = ',istop
    print *, '(REPORTED BY SOLVER) ierr                             = ',ierr
    print *, '(REPORTED BY SOLVER) f                                = ',f
    print *, '(REPORTED BY SOLVER) csupn                            = ',csupn
    print *, '(REPORTED BY SOLVER) ssupn                            = ',ssupn
    print *, '(REPORTED BY SOLVER) nlpsupn                          = ',nlpsupn
    print *, '(REPORTED BY SOLVER) bounds violation                 = ',bdsvio
    print *, '(REPORTED BY SOLVER) Number of outer iterations       = ',outiter
    print *, '(REPORTED BY SOLVER) Number of inner iterations       = ',totiter
    print *, '(REPORTED BY SOLVER) Number of Newton-KKT trials      = ',nwcalls
    print *, '(REPORTED BY SOLVER) Number of Newton-KKT iterations  = ',nwtotit

    print *, ''
    print *, '(COMPUTED BY CALLER) Number of calls to evalf         = ',pdata%counters(1)
    print *, '(COMPUTED BY CALLER) Number of calls to evalg         = ',pdata%counters(2)
    print *, '(COMPUTED BY CALLER) Number of calls to evalc         = ',pdata%counters(3)
    print *, '(COMPUTED BY CALLER) Number of calls to evalj         = ',pdata%counters(4)
    print *, '(COMPUTED BY CALLER) Number of calls to evalhl        = ',pdata%counters(5)
    print *, '(COMPUTED BY CALLER) CPU time in seconds              = ',finish - start

    ! *****************************************************************
    ! *****************************************************************
    ! Just checking ...

    inform = 0

    ! call user_evalf(n,x,f,inform,c_loc(pdata))

    if ( inform .ne. 0 ) then
      print *, 'error when calling evalf in the main file. '
      stop
    end if

    ! call user_evalc(n,x,m,p,c,inform,c_loc(pdata))

    if ( inform .ne. 0 ) then
      print *, 'error when calling evalc in the main file. '
      stop
    end if

    csupn = max( 0.0d0, max( maxval( abs( c(1:m) ) ), maxval( c(m+1:m+p) ) ) )

    bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

    print *, ''
    print *, '(COMPUTED BY CALLER) f                                = ',f
    print *, '(COMPUTED BY CALLER) csupn                            = ',csupn
    print *, '(COMPUTED BY CALLER) bounds violation                 = ',bdsvio

    print *, ''
    print *, 'When a quantity appears as computed by solver and computed by caller, they must coincide.'
    print *, '(In case they do not coincide, please report it as a bug.)'
    ! *****************************************************************
    ! *****************************************************************

    write (*,*) "     x = ", x
    write (*,*) "lambda = ", lambda

    deallocate(c,stat=allocerr)

    if ( allocerr .ne. 0 ) then
      print *, 'Deallocation error.'
      stop
    end if

    stop
  end subroutine init

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf_original(n,x,f,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! This routine must compute the objective function.

    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1

    f = ( x(1) + 3.0d0 * x(2) + x(3) ) ** 2.0d0 + 4.0d0 * (x(1) - x(2)) ** 2.0d0

  end subroutine evalf_original

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg_original(n,x,g,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! This routine must compute the gradient of the objective
    ! function.

    ! LOCAL SCALARS
    real(kind=8) :: t1, t2
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1

    t1 = x(1) + 3.0d0 * x(2) + x(3)
    t2 = x(1) - x(2)
    g(1) = 2.0d0 * t1 + 8.0d0 * t2
    g(2) = 6.0d0 * t1 - 8.0d0 * t2
    g(3) = 2.0d0 * t1

  end subroutine evalg_original

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalc_original(n,x,m,p,c,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m+p)

    ! This routine must compute all the m+p constraints.

    ! LOCAL SCALARS
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    c(1) = 1.0d0 - x(1) - x(2) - x(3)
    c(2) = - 6.0d0 * x(2) - 4.0d0 * x(3) + x(1) ** 3.0d0 + 3.0d0

  end subroutine evalc_original

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalj_original(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    logical, intent(in) :: ind(m+p)
    real(kind=8), intent(in) :: x(n)
    logical, intent(out) :: sorted(m+p)
    integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
    real(kind=8), intent(out) :: jval(lim)

    ! This routine must compute the Jacobian of the constraints. In
    ! fact, only gradients of constraints j such that ind(j) =
    ! .true. need to be computed.

    ! LOCAL SCALARS
    integer :: i, nnz1, nnz2
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(4) = pdata%counters(4) + 1

    ! Only gradients of constraints j such that ind(j) = .true. need
    ! to be computed.

    nnz1 = 0
    nnz2 = 0

    if ( ind(1) ) then
       if ( lim .lt. n ) then
          inform = -94
          return
       end if

       jsta(1) = 1
       jlen(1) = n
       nnz1 = nnz1 + 1

       jvar(1:n) = (/ (i,i=1,n) /)
       jval(1:n) = -1.0d0
       nnz2 = nnz2 + n

       ! Says whether the variables' indices in jvar (related to this
       ! constraint) are in increasing order. In case they are,
       ! Algencan takes advantage of this. Implement sorted gradients
       ! of constraints if you can do this in a natural (cheap)
       ! way. Under no circumnstance use a sorting algorithm. (It is
       ! better to set sorted(1) = .false. in this case.)

       sorted(1) = .true.
    end if

    if ( ind(2) ) then
       if ( lim .lt. n ) then
          inform = -94
          return
       end if

       jsta(nnz1+1) = nnz2 + 1
       jlen(nnz1+1) = n

       jvar(nnz2+1:nnz2+n) = (/ (i,i=1,n) /)

       jval(nnz2+1) = - 3.0d0 * x(1) ** 2.0d0
       jval(nnz2+2) = 6.0d0
       jval(nnz2+3) = 4.0d0

       sorted(2) = .true.
    end if

  end subroutine evalj_original

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhl_original(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in) :: inclf
    integer, intent(in) :: m,n,lim,p
    integer, intent(out) :: hlnnz
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m+p),x(n)
    integer, intent(out) :: hlrow(lim),hlcol(lim)
    real(kind=8), intent(out) :: hlval(lim)

    ! This routine must compute the Hessian of the Lagrangian. The
    ! Hessian of the objective function must NOT be included if inclf
    ! = .false.

    ! LOCAL SCALARS
    integer :: i
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(5) = pdata%counters(5) + 1

    hlnnz = 0

    ! If .not. inclf then the Hessian of the objective function must not be included

    if ( inclf ) then
       if ( hlnnz + 2 .gt. lim ) then
          inform = -95
          return
       end if

       hlnnz = 6

       hlrow(1:6) = [      1,      2,      2,     3,     3,     3 ]
       hlcol(1:6) = [      1,      1,      2,     1,     2,     3 ]
       hlval(1:6) = [ 10.0d0, -2.0d0, 26.0d0, 2.0d0, 6.0d0, 2.0d0 ]

    end if

    ! Note that entries of the Hessian of the Lagrangian can be
    ! repeated. If this is case, them sum of repeated entrances is
    ! considered. This feature simplifies the construction of the
    ! Hessian of the Lagrangian.

    if ( hlnnz + 1 .gt. lim ) then
       inform = -95
       return
    end if

    hlnnz = hlnnz + 1

    hlrow(hlnnz) = 1
    hlcol(hlnnz) = 1
    hlval(hlnnz) = lambda(2) * ( - 6.0d0 * x(1) )

  end subroutine evalhl_original

end module algencanma
