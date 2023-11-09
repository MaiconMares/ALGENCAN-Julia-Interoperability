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

      implicit none

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
    end subroutine evalg

    subroutine evalc(n,x,m,p,c,inform,pdataptr)
      import c_ptr

      implicit none

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

    end subroutine evalj

    subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)
      import c_ptr
      
      implicit none

      logical, intent(in) :: inclf
      integer, intent(in) :: m,n,lim,p
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
    x,n,f,lind,lbnd,uind,ubnd,m,p,lambda,jnnzmax,hlnnzmax,epsfeas,epscompl,  &
    epsopt,rhoauto,rhoini,scale,extallowed,corrin,inform,maxoutit,pdata)

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(evalf) :: user_evalf
    procedure(evalg) :: user_evalg
    procedure(evalc) :: user_evalc
    procedure(evalj) :: user_evalj
    procedure(evalhl) :: user_evalhl

    ! PARAMETERS
    integer, target, intent(inout) :: n, m, p
    logical, intent(in) :: corrin,extallowed,rhoauto,scale, lind(n),uind(n)
    real(kind=8), intent(inout) :: x(n), epsfeas,epscompl,epsopt,rhoini
    real(kind=8), target, intent(inout) :: lambda(m+p)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    integer, intent(in) :: hlnnzmax, jnnzmax
    real(kind=8), target, intent(out) :: f
    integer, target, intent(inout) :: inform
    integer, intent(in) :: maxoutit
    type(pdata_type), target, intent(inout) :: pdata

    
    ! LOCAL SCALARS
    integer :: allocerr,ierr,istop,nwcalls,nwtotit,outiter,totiter
    real(kind=8) :: bdsvio,csupn,finish,nlpsupn,ssupn,start
    !type(pdata_type), target :: pdata
    real(kind=8) :: c(m+p)
    
    if ( allocerr .ne. 0 ) then
      print *, 'Allocation error.'
      stop
    end if
    
    if ( allocerr .ne. 0 ) then
      print *, 'Allocation error.'
      stop
    end if

    call cpu_time(start)
    
    call algencan(user_evalf,user_evalg,user_evalc,user_evalj,user_evalhl,jnnzmax,hlnnzmax, &
        n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
        scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
        outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

    call cpu_time(finish)

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

    call user_evalf(n,x,f,inform,c_loc(pdata))

    if ( inform .ne. 0 ) then
      print *, 'error when calling evalf in the main file. '
      stop
    end if

    call user_evalc(n,x,m,p,c,inform,c_loc(pdata))

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

    write (*,*) "lambda = ", lambda

    if ( allocerr .ne. 0 ) then
      print *, 'Deallocation error.'
      stop
    end if

    !stop
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

    f = ( x(1) - 5.0d0 ) ** 2.0d0 + x(2)**2 - 25.0d0

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
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1

    g(1) = 2*x(1) - 10.0d0
    g(2) = 2*x(2)

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

    c(1) = x(1)**2 - x(2)

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
       jval(1) = 2.0d0 * x(1)
       jval(2) = -1.0d0
       nnz2 = nnz2 + n

       ! Says whether the variables' indices in jvar (related to this
       ! constraint) are in increasing order. In case they are,
       ! Algencan takes advantage of this. Implement sorted gradients
       ! of constraints if you can do this in a natural (cheap)
       ! way. Under no circumnstance use a sorting algorithm. (It is
       ! better to set sorted(1) = .false. in this case.)

       sorted(1) = .true.
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

       hlnnz = 3

       hlrow(1:3) = [      1,      2,      1]
       hlcol(1:3) = [      1,      2,      1]
       hlval(1:3) = [ 2.0d0, 2.0d0, 2.0d0*lambda(1)]

      print *, "hlrow = ", hlrow(1:10)
      print *, "hlcol = ", hlcol(1:10)
      print *, "hlval = ", hlval(1:10)
      print *, "lambda = ", lambda

    end if

    ! Note that entries of the Hessian of the Lagrangian can be
    ! repeated. If this is case, them sum of repeated entrances is
    ! considered. This feature simplifies the construction of the
    ! Hessian of the Lagrangian.

    if ( hlnnz + 1 .gt. lim ) then
       inform = -95
       return
    end if

  end subroutine evalhl_original

end module algencanma
