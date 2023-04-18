! *****************************************************************
! *****************************************************************

program algencama

  use bmalgencan, only: algencan
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer

  implicit none

  type :: pdata_type
     integer :: counters(5) = 0

     integer :: m_orig,n_orig,nfree
     logical, allocatable :: free(:),frelem(:)
     integer, allocatable :: ccor(:),cmap(:),freeind(:),ifrelem(:),jperm(:),renumbered(:)
     real(kind=8), allocatable :: ca(:),cb(:),gfull(:),v(:),xfull(:)
  end type pdata_type

  ! LOCAL SCALARS
  logical :: corrin,extallowed,rhoauto,scale
  integer :: allocerr,e_order,hlnnzmax,i,ierr,istop,j,jnnz,jnnzmax,jfun(1),jvar(1),l_order, &
       v_order,m,maxoutit,n,ncdouble,nwcalls,nwtotit,outiter,p,status,totiter
  real(kind=8) :: bdsvio,csupn,epsfeas,epscompl,epsopt,f,finish,jval(1),nlpsupn,rhoini,ssupn,start
  type(pdata_type), target :: pdata

  ! LOCAL ARRAYS
  character(len=10) :: pname
  logical, allocatable :: cdouble(:),equatn(:),linear(:),clind(:),cuind(:),lind(:),uind(:)
  real(kind=8), allocatable :: c(:),clbnd(:),cubnd(:),lbnd(:),lbnd_orig(:),ubnd(:),ubnd_orig(:), &
       lambda(:),lambda_orig(:),x(:)

  open(10,file='OUTSDIF.d',form='formatted',status='old')

  rewind 10

  call cutest_cdimen(status,10,pdata%n_orig,pdata%m_orig)

  if ( status .ne. 0 ) then
     write(*,*) 'cutest_cdimen: status not equal to zero. status = ',status
     stop
  end if

  allocate(lbnd_orig(pdata%n_orig),ubnd_orig(pdata%n_orig),pdata%xfull(pdata%n_orig), &
       pdata%gfull(pdata%n_orig),pdata%free(pdata%n_orig),pdata%renumbered(pdata%n_orig), &
       lambda_orig(pdata%m_orig),clind(pdata%m_orig),clbnd(pdata%m_orig),cuind(pdata%m_orig), &
       cubnd(pdata%m_orig),equatn(pdata%m_orig),linear(pdata%m_orig),cdouble(pdata%m_orig), &
       pdata%v(pdata%m_orig),pdata%ccor(pdata%m_orig),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  e_order = 1
  l_order = 0
  v_order = 0

  call CUTEST_csetup(status,10,6,20,pdata%n_orig,pdata%m_orig,pdata%xfull,lbnd_orig,ubnd_orig, &
       lambda_orig,clbnd,cubnd,equatn,linear,e_order,l_order,v_order)

  if ( status .ne. 0 ) then
     write(*,*) 'cutest_csetup: status not equal to zero. status = ',status
     stop
  end if

  where ( lbnd_orig(1:pdata%n_orig) .eq. ubnd_orig(1:pdata%n_orig) )
     pdata%free(1:pdata%n_orig) = .false.
     pdata%xfull(1:pdata%n_orig) = lbnd_orig(1:pdata%n_orig)
  elsewhere
     pdata%free(1:pdata%n_orig) = .true.
  end where

  pdata%nfree = count( pdata%free(1:pdata%n_orig) )

  allocate(lind(pdata%nfree),lbnd(pdata%nfree),uind(pdata%nfree),ubnd(pdata%nfree),x(pdata%nfree), &
       pdata%freeind(pdata%nfree),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.',allocerr
     stop
  end if

  pdata%freeind(1:pdata%nfree) = pack( (/ (i,i=1,pdata%n_orig) /), pdata%free(1:pdata%n_orig) )
  pdata%renumbered(1:pdata%n_orig) = unpack( (/ (i,i=1,pdata%nfree) /), pdata%free(1:pdata%n_orig), 0 )

  n = pdata%nfree

  x(1:n) = pdata%xfull(pdata%freeind(1:pdata%nfree))
  lbnd(1:n) = lbnd_orig(pdata%freeind(1:pdata%nfree))
  ubnd(1:n) = ubnd_orig(pdata%freeind(1:pdata%nfree))

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

  deallocate(lbnd_orig,ubnd_orig,stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Dellocation error.'
     stop
  end if

  m = count( equatn(1:pdata%m_orig) )

  if ( any( .not. equatn(1:m) ) .or. any( equatn(m+1:pdata%m_orig) ) ) then
     write(*,*) 'There is something wrong with equatn, cl, and cu in the problem definition.'
     stop
  end if

  if ( any( clbnd(1:m) .ne. cubnd(1:m) ) ) then
     write(*,*) 'There is something wrong with equatn, cl, and cu in the problem definition.'
     stop
  end if

  where( clbnd(m+1:pdata%m_orig) .gt. - 1.0d+20 )
     clind(m+1:pdata%m_orig) = .true.
  elsewhere
     clind(m+1:pdata%m_orig) = .false.
  end where

  where( cubnd(m+1:pdata%m_orig) .lt.   1.0d+20 )
     cuind(m+1:pdata%m_orig) = .true.
  elsewhere
     cuind(m+1:pdata%m_orig) = .false.
  end where

  if ( any( .not. clind(m+1:pdata%m_orig) .and. .not. cuind(m+1:pdata%m_orig) ) ) then
     write(*,*) 'There is something wrong with equatn, cl, and cu in the problem definition.'
     stop
  end if

  cdouble(m+1:pdata%m_orig) = clind(m+1:pdata%m_orig) .and. cuind(m+1:pdata%m_orig)

  ncdouble = count( cdouble(m+1:pdata%m_orig) )

  p = pdata%m_orig - m + ncdouble

  allocate(lambda(m+p),c(m+p),pdata%ca(m+p),pdata%cb(m+p),pdata%cmap(m+p),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  pdata%cmap(1:pdata%m_orig) = (/ (j,j=1,pdata%m_orig) /)
  pdata%cmap(pdata%m_orig+1:m+p) = pack( (/ (j,j=m+1,pdata%m_orig) /), cdouble(m+1:pdata%m_orig) )

  pdata%ccor(1:pdata%m_orig) = 0
  pdata%ccor(pdata%cmap(pdata%m_orig+1:m+p)) = (/ (j,j=pdata%m_orig+1,m+p) /)

  pdata%ca(1:m) = 1.0d0
  pdata%cb(1:m) = -cubnd(1:m)

  where( cuind(m+1:pdata%m_orig) )
     pdata%ca(m+1:pdata%m_orig)   = 1.0d0
     pdata%cb(m+1:pdata%m_orig)   = -cubnd(m+1:pdata%m_orig)
  elsewhere
     pdata%ca(m+1:pdata%m_orig)   = -1.0d0
     pdata%cb(m+1:pdata%m_orig)   = clbnd(m+1:pdata%m_orig)
  end where

  where( cdouble(m+1:pdata%m_orig) )
     pdata%ca(pdata%m_orig+1:m+p) = -1.0d0
     pdata%cb(pdata%m_orig+1:m+p) = clbnd(pdata%cmap(pdata%m_orig+1:m+p))
  end where

  deallocate(lambda_orig,clind,clbnd,cuind,cubnd,equatn,linear,cdouble,stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if

  call cutest_cdimsj(status,jnnzmax)

  if ( status .ne. 0 ) then
     write(*,*) 'cutest_cdimsj: status not equal to zero. status = ',status
     stop
  end if

  jnnzmax = 2 * jnnzmax

  hlnnzmax = huge( 1 )

  allocate(pdata%jperm(jnnzmax),pdata%frelem(max(jnnzmax,hlnnzmax)),pdata%ifrelem(max(jnnzmax,hlnnzmax)),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error.'
     stop
  end if

  call cutest_pname(status,10,pname)

  if ( status .ne. 0 ) then
     write(*,*) 'cutest_pname: status not equal to zero. status = ',status
     stop
  end if

  close(10)

  lambda(1:m+p) = 0.0d0

  epsfeas  = 1.0d-08
  epscompl = 1.0d-08
  epsopt   = 1.0d-08

  maxoutit = 50

  rhoauto = .true.

  if ( .not. rhoauto ) then
     rhoini = 1.0d-08
  end if

  scale = .false.

  extallowed = .true.

  corrin = .false.

  call cpu_time(start)

  call algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
       n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
       scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn,bdsvio, &
       outiter,totiter,nwcalls,nwtotit,ierr,istop,c_loc(pdata))

  call cpu_time(finish)

  write(*,*)
  write(*,*) 'Problem name                                          = ',pname
  write(*,*) 'Number of variables                                   = ',n
  write(*,*) 'Number of equality constraints                        = ',m
  write(*,*) 'Number of one-sided inequalities                      = ',pdata%m_orig - m - ncdouble
  write(*,*) 'Number of two-sided inequalities                      = ',ncdouble

  write(*,*)
  write(*,*) '(REPORTED BY SOLVER) istop                            = ',istop
  write(*,*) '(REPORTED BY SOLVER) ierr                             = ',ierr
  write(*,*) '(REPORTED BY SOLVER) f                                = ',f
  write(*,*) '(REPORTED BY SOLVER) csupn                            = ',csupn
  write(*,*) '(REPORTED BY SOLVER) ssupn                            = ',ssupn
  write(*,*) '(REPORTED BY SOLVER) nlpsupn                          = ',nlpsupn
  write(*,*) '(REPORTED BY SOLVER) bounds violation                 = ',bdsvio
  write(*,*) '(REPORTED BY SOLVER) Number of outer iterations       = ',outiter
  write(*,*) '(REPORTED BY SOLVER) Number of inner iterations       = ',totiter
  write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT trials      = ',nwcalls
  write(*,*) '(REPORTED BY SOLVER) Number of Newton-KKT iterations  = ',nwtotit

  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalf         = ',pdata%counters(1)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalg         = ',pdata%counters(2)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalc         = ',pdata%counters(3)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalj         = ',pdata%counters(4)
  write(*,*) '(COMPUTED BY CALLER) Number of calls to evalhl        = ',pdata%counters(5)
  write(*,*) '(COMPUTED BY CALLER) CPU time in seconds              = ',finish - start

  ! *****************************************************************
  ! *****************************************************************
  ! Just checking ...

  pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

  call cutest_cofg(status,pdata%n_orig,pdata%xfull,f,pdata%gfull,.false.)

  if ( status .ne. 0 ) then
     write(*,*) 'error when calling cutest_cofg in the main file. '
     stop
  end if

  call cutest_ccfsg(status,pdata%n_orig,pdata%m_orig,pdata%xfull,c,jnnz,1,jval,jvar,jfun,.false.)

  if ( status .ne. 0 ) then
     write(*,*) 'error when calling cutest_ccfsg in the main file. '
     stop
  end if

  c(1:m+p) = pdata%ca(1:m+p) * c(pdata%cmap(1:m+p)) + pdata%cb(1:m+p)

  csupn = max( 0.0d0, max( maxval( abs( c(1:m) ) ), maxval( c(m+1:m+p) ) ) )

  bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

  write(*,*)
  write(*,*) '(COMPUTED BY CALLER) f                                = ',f
  write(*,*) '(COMPUTED BY CALLER) csupn                            = ',csupn
  write(*,*) '(COMPUTED BY CALLER) bounds violation                 = ',bdsvio

  write(*,*)
  write(*,*) 'When a quantity appears as computed by solver and computed by caller, they must coincide.'
  write(*,*) '(In case they do not coincide, please report it as a bug.)'
  ! *****************************************************************
  ! *****************************************************************

  open(20,file='tabline.txt')
  write(20,9000) istop,ierr,finish - start,n,m+p,f,csupn,ssupn,nlpsupn,bdsvio, &
       outiter,totiter,nwcalls,nwtotit,pdata%counters(1:5)
  close(20)

  deallocate(lind,lbnd,uind,ubnd,x,lambda,c,pdata%v,pdata%ca,pdata%cb,pdata%cmap, &
       pdata%ccor,pdata%jperm,pdata%frelem,pdata%ifrelem,pdata%free,pdata%freeind, &
       pdata%renumbered,pdata%xfull,pdata%gfull,stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error.'
     stop
  end if

  stop

9000 format(1X,I2,1X,I3,0P,F12.6,2(1X,I6),1X,1P,D24.16,4(1X,1P,D7.1),9(1X,I8))

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

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(1) = pdata%counters(1) + 1

    pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

    call cutest_cofg(status,pdata%n_orig,pdata%xfull,f,pdata%gfull,.false.)
    if ( status .ne. 0 ) then
       inform = -91
       return
    end if

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
    real(kind=8) :: f
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(2) = pdata%counters(2) + 1

    pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

    call cutest_cofg(status,pdata%n_orig,pdata%xfull,f,pdata%gfull,.true.)

    where ( isnan( pdata%gfull(1:pdata%n_orig) ) )
       pdata%gfull(1:pdata%n_orig) = 0.0d0
    end where

    ! Remove elements related to fixed variables

    g(1:n) = pdata%gfull(pdata%freeind(1:pdata%nfree))

    if ( status .ne. 0 ) then
       inform = -92
       return
    end if

  end subroutine evalg

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalc(n,x,m,p,c,inform,pdataptr)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: inform
    type(c_ptr), optional, intent(in) :: pdataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: c(m+p)

    ! LOCAL SCALARS
    integer :: jfun(1),jvar(1),jnnz,status
    real(kind=8) :: jval(1)
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(3) = pdata%counters(3) + 1

    pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

    call cutest_ccfsg(status,pdata%n_orig,pdata%m_orig,pdata%xfull,c,jnnz,1,jval,jvar,jfun,.false.)
    if ( status .ne. 0 ) then
       inform = -93
       return
    end if

    c(1:m+p) = pdata%ca(1:m+p) * c(pdata%cmap(1:m+p)) + pdata%cb(1:m+p)

  end subroutine evalc

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalj(n,x,m,p,ind,sorted,jsta,jlen,lim,jvar,jval,inform,pdataptr)

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

    ! LOCAL SCALARS
    integer :: jnnz,nfrelem,nnew
    type(pdata_type), pointer :: pdata

    ! LOCAL ARRAYS
    logical :: tmp(lim)
    real(kind=8) :: c(m+p) ! m+p is an upper bound of m_orig
    integer :: jfun(lim)

    ! Only the gradients of the constraints indicated in the logical
    ! array ind should be computed. This fact is being ignored here
    ! and all of them are being computed.

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(4) = pdata%counters(4) + 1

    pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

    call CUTEST_ccfsg(status,pdata%n_orig,pdata%m_orig,pdata%xfull,c,jnnz,lim,jval,jvar,jfun,.true.)
    if ( status .ne. 0 ) then
       inform = -94
       return
    end if

    where ( isnan( jval(1:jnnz) ) )
       jval(1:jnnz) = 0.0d0
    end where

    ! Remove elements related to fixed variables

    if ( pdata%nfree .ne. pdata%n_orig ) then
       pdata%frelem(1:jnnz) = pdata%free(jvar(1:jnnz))

       nfrelem = count( pdata%frelem(1:jnnz) )
       pdata%ifrelem(1:nfrelem) = pack( (/ (i, i=1,jnnz) /), pdata%frelem(1:jnnz) )

       jnnz = nfrelem
       jfun(1:nfrelem) = jfun(pdata%ifrelem(1:nfrelem))
       jvar(1:nfrelem) = pdata%renumbered(jvar(pdata%ifrelem(1:nfrelem)))
       jval(1:nfrelem) = jval(pdata%ifrelem(1:nfrelem))
    end if

    ! Duplicate two-sides constraints

    tmp(1:jnnz) = pdata%ccor(jfun(1:jnnz)) .ne. 0

    nnew = count( tmp(1:jnnz) )

    if ( jnnz + nnew .gt. lim ) then
       inform = -96
       return
    end if

    jfun(jnnz+1:jnnz+nnew) = pack( pdata%ccor(jfun(1:jnnz)), tmp(1:jnnz) )
    jvar(jnnz+1:jnnz+nnew) = pack( jvar(1:jnnz), tmp(1:jnnz) )
    jval(jnnz+1:jnnz+nnew) = pack( jval(1:jnnz), tmp(1:jnnz) )

    jnnz = jnnz + nnew

    jval(1:jnnz) = pdata%ca(jfun(1:jnnz)) * jval(1:jnnz)

    ! Convert to compressed sparse row format

    call coo2csr(m+p,jnnz,jfun,jvar,jval,jlen,jsta)

    sorted(1:m+p) = .false.

  end subroutine evalj

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhl(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,inform,pdataptr)

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

    ! LOCAL SCALARS
    integer :: nfrelem,status
    type(pdata_type), pointer :: pdata

    call c_f_pointer(pdataptr,pdata)
    pdata%counters(5) = pdata%counters(5) + 1

    pdata%xfull(pdata%freeind(1:pdata%nfree)) = x(1:n)

    pdata%v(1:pdata%m_orig) = 0.0d0
    pdata%v(pdata%cmap(1:m+p)) = pdata%v(pdata%cmap(1:m+p)) + pdata%ca(1:m+p) * lambda(1:m+p)

    if ( inclf ) then
       call cutest_csh(status,pdata%n_orig,pdata%m_orig,pdata%xfull,pdata%v,hlnnz,lim,hlval,hlrow,hlcol)
    else
       call cutest_cshc(status,pdata%n_orig,pdata%m_orig,pdata%xfull,pdata%v,hlnnz,lim,hlval,hlrow,hlcol)
    end if

    where ( isnan( hlval(1:hlnnz) ) )
       hlval(1:hlnnz) = 0.0d0
    end where

    if ( status .ne. 0 ) then
       inform = -95
       return
    end if

    ! Remove elements related to fixed variables

    if ( pdata%nfree .ne. pdata%n_orig ) then
       pdata%frelem(1:hlnnz) = pdata%free(hlrow(1:hlnnz)) .and. pdata%free(hlcol(1:hlnnz))

       nfrelem = count( pdata%frelem(1:hlnnz) )
       pdata%ifrelem(1:nfrelem) = pack( (/ (i, i=1,hlnnz) /), pdata%frelem(1:hlnnz) )

       hlnnz = nfrelem
       hlrow(1:nfrelem) = pdata%renumbered(hlrow(pdata%ifrelem(1:nfrelem)))
       hlcol(1:nfrelem) = pdata%renumbered(hlcol(pdata%ifrelem(1:nfrelem)))
       hlval(1:nfrelem) = hlval(pdata%ifrelem(1:nfrelem))
    end if

  end subroutine evalhl

! ******************************************************************
! ******************************************************************

  subroutine coo2csr(nrow,nnz,arow,acol,aval,alen,asta)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nrow,nnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: acol(nnz),alen(nrow),arow(nnz),asta(nrow)
    real(kind=8), intent(inout) :: aval(nnz)

    ! This subroutines converts a matrix from coordinate format to
    ! compressed sparse row format.

    ! LOCAL SCALARS
    integer :: i,j,col,coltmp,row,rowtmp
    real(kind=8) :: val,valtmp

    alen(1:nrow) = 0

    do i = 1,nnz
       row = arow(i)
       alen(row) = alen(row) + 1
    end do

    asta(1) = 1
    do i = 2,nrow
       asta(i) = asta(i-1) + alen(i-1)
    end do

    do i = 1,nnz
       val = aval(i)
       col = acol(i)
       row = arow(i)

       arow(i) = - 1

       do while ( row .ge. 0 )
          j = asta(row)
          asta(row) = j + 1

          valtmp = aval(j)
          coltmp = acol(j)
          rowtmp = arow(j)

          aval(j) = val
          acol(j) = col
          arow(j) = - 1

          val = valtmp
          col = coltmp
          row = rowtmp
       end do
    end do

    asta(1:nrow) = asta(1:nrow) - alen(1:nrow)

  end subroutine coo2csr

end program algencama

