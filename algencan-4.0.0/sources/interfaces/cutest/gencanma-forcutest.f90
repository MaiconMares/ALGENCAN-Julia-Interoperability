! *****************************************************************
! *****************************************************************

program gencanma

  use bmgencan, only: gencan, genunc
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

  implicit none

  type pdata_type
     integer :: counters(3) = 0
  end type pdata_type

  ! LOCAL SCALARS
  logical :: extallowed,hfixstr
  integer :: allocerr,hnnzmax,ierr,istop,iter,maxit,n,nbds,status
  real(kind=8) :: bdsvio,eps,f,finish,ftarget,gpsupn,start
  type(pdata_type), target :: pdata

  ! LOCAL ARRAYS
  character(len=10) :: pname
  logical, allocatable :: lind(:),uind(:)
  real(kind=8), allocatable :: g(:),lbnd(:),ubnd(:),x(:)

  open(10,file='OUTSDIF.d',form='formatted',status='old')

  rewind 10

  call cutest_udimen(status,10,n)
  if ( status .ne. 0 ) then
     write(*,*) 'cutest_udimen: status not equal to zero. status = ',status
     stop
  end if

  allocate(g(n),lind(n),lbnd(n),uind(n),ubnd(n),x(n),stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  call cutest_usetup(status,10,6,4,n,x,lbnd,ubnd)
  if ( status .ne. 0 ) then
     write(*,*) 'cutest_usetup: status not equal to zero. status = ',status
     stop
  end if

  where( lbnd(1:n) .gt. - 1.0d+20 )
     lind(1:n) = .true.
  elsewhere
     lind(1:n) = .false.
  end where

  where( ubnd(1:n) .lt.   1.0d+20 )
     uind(1:n) = .true.
  elsewhere
     uind(1:n) = .false.
  end where

  call cutest_pname(status,10,pname)

  call cutest_udimsh(status,hnnzmax)
  if ( status .ne. 0 ) then
     write(*,*) 'cutest_udimsh: status not equal to zero. status = ',status
     stop
  end if

  hnnzmax = hnnzmax + n
  
  close(10)

  ftarget     =  - 1.0d+12
  eps         =    1.0d-08
  maxit       =      50000
  hfixstr     =     .true.
  extallowed  =     .true.

  nbds = count( lind(1:n) ) + count( uind(1:n) )

  call cpu_time(start)

  if ( nbds .eq. 0 ) then
     write(*,*) 'The problem is unconstrained.'
     call genunc(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,f,g,gpsupn, &
          ftarget,eps,maxit,extallowed,iter,ierr,istop,pdata=c_loc(pdata))
     
  else
     write(*,*) 'The problem is a bound constrained problem.'
     call gencan(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,lind,lbnd, &
          uind,ubnd,f,g,gpsupn,ftarget,eps,maxit,extallowed,iter,ierr, &
          istop,pdata=c_loc(pdata))
  end if
  
  call cpu_time(finish)

  write(*,*)
  write(*,*) 'Problem name                                          = ',pname
  write(*,*) 'Number of variables                                   = ',n
  write(*,*) 'Number of lower and upper bounds                      = ',count( lind(1:n) ) + count( uind(1:n) )
  write(*,*)
  write(*,*) '(REPORTED BY SOLVER) istop                            = ',istop
  write(*,*) '(REPORTED BY SOLVER) ierr                             = ',ierr
  write(*,*) '(REPORTED BY SOLVER) f                                = ',f
  write(*,*) '(REPORTED BY SOLVER) gpsupn                           = ',gpsupn
  write(*,*) '(REPORTED BY SOLVER) Number of iterations             = ',iter
  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) Number of function evaluations   = ',pdata%counters(1)
  write(*,*) '(COMPUTED BY CALLER) Number of gradient evaluations   = ',pdata%counters(2)
  write(*,*) '(COMPUTED BY CALLER) Number of Hessian  evaluations   = ',pdata%counters(3)
  write(*,*) '(COMPUTED BY CALLER) CPU time in seconds              = ',finish - start

  ! *****************************************************************
  ! *****************************************************************
  ! Just checking ...
  call cutest_ufn(status,n,x,f)
  if ( status .ne. 0 ) then
     write(*,*) 'error when calling cutest_ufn in the main file.'
     stop
  end if
  
  bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) f                                = ',f
  write(*,*) '(COMPUTED BY CALLER) bounds violation                 = ',bdsvio

  write(*,*)
  write(*,*) 'When a quantity appears as computed by solver and computed by caller, they must coincide.'
  write(*,*) '(In case they do not coincide, please report it as a bug.)'
  ! *****************************************************************
  ! *****************************************************************
  
  open(20,file='tabline.txt')
  write(20,9000) istop,ierr,finish - start,n,f,gpsupn,iter,pdata%counters(1:3)
  close(20)

  deallocate(g,lind,lbnd,uind,ubnd,x,stat=allocerr)
  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if

  stop

9000 format(1X,I2,1X,I3,0P,F12.6,1X,I6,1X,1P,D24.16,1X,1P,D7.1,4(1X,I8))

contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalf(n,x,f,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    real(kind=8), intent(out) :: f
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    
    ! LOCAL SCALARS
    integer :: status
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr, pdata)
    pdata%counters(1) = pdata%counters(1) + 1
    
    call cutest_ufn(status,n,x,f)
    if ( status .ne. 0 ) inform = -91
  
  end subroutine evalf

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalg(n,x,g,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr
  
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: g(n)
    
    ! LOCAL SCALARS
    integer :: status
    type(pdata_type), pointer :: pdata
    
    call c_f_pointer(pdataptr, pdata)
    pdata%counters(2) = pdata%counters(2) + 1
    
    call cutest_ugr(status,n,x,g)
    if ( status .ne. 0 ) inform = -92

    where ( isnan( g(1:n) ) )
       g(1:n) = 0.0d0
    end where
    
  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalh(n,x,lim,hnnz,hrow,hcol,hval,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,n
    integer, intent(inout) :: inform
    integer, intent(out) :: hnnz
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    integer, intent(out) :: hcol(lim),hrow(lim)
    real(kind=8), intent(out) :: hval(lim)

    ! LOCAL SCALARS
    integer :: status
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr, pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    ! Upper triangle
    call cutest_ush(status,n,x,hnnz,lim,hval,hrow,hcol)
    if ( status .ne. 0 ) inform = -93
    
    where ( isnan( hval(1:hnnz) ) )
       hval(1:hnnz) = 0.0d0
    end where
    
  end subroutine evalh

end program gencanma

