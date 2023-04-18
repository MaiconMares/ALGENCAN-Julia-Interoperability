module memev

  use iso_c_binding, only: c_ptr
  
  implicit none

  interface
     subroutine fsub(n,x,f,ierr,pdata)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: f 
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine fsub

     subroutine gsub(n,x,g,ierr,pdata)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: g(n)
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine gsub

     subroutine csub(n,x,m,p,c,ierr,pdata)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: m,n,p
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: c(m+p) 
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine csub

     subroutine jsub(n,x,m,p,ind,jsorted,jsta,jlen,lim,jvar,jval,ierr,pdata)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: lim,m,n,p
       logical, intent(in) :: ind(m+p)
       real(kind=8), intent(in) :: x(n)
       logical, intent(out) :: jsorted(m+p)
       integer, intent(out) :: jsta(m+p),jlen(m+p),jvar(lim)
       real(kind=8), intent(out) :: jval(lim)
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine jsub

     subroutine hlsub(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,ierr,pdata)
       use iso_c_binding, only: c_ptr
       implicit none
       logical, intent(in) :: inclf
       integer, intent(in) :: m,n,lim,p
       real(kind=8), intent(in) :: lambda(m+p),x(n)
       integer, intent(out) :: hlnnz
       integer, intent(out) :: hlrow(lim),hlcol(lim)
       real(kind=8), intent(out) :: hlval(lim)
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine hlsub
  end interface

  type :: memev_type
     procedure(fsub), pointer, nopass :: evalf
     procedure(gsub), pointer, nopass :: evalg
     procedure(csub), pointer, nopass :: evalc
     procedure(jsub), pointer, nopass :: evalj
     procedure(hlsub), pointer, nopass :: evalhl
     
     logical :: pdatapresent,scale
     integer :: jnnzmax,m,n,p
     real(kind=8) :: sf
     type(c_ptr), pointer :: pdata
     
     real(kind=8), allocatable :: sc(:)
     
     logical :: firstf,firstg,firstc,firstj
     real(kind=8) :: f,scaledf
     
     logical, allocatable :: jind(:),jsorted(:)
     integer, allocatable :: jsta(:),jlen(:),jvar(:)
     real(kind=8), allocatable :: c(:),g(:),jval(:),scaledc(:),scaledg(:), &
          scaledjval(:),xf(:),xg(:),xc(:),xj(:)
  end type memev_type
  
  type :: memev_backup_type
     logical :: firstf,firstg,firstc,firstj
     real(kind=8) :: f,scaledf
     
     logical, allocatable :: jind(:),jsorted(:)
     integer, allocatable :: jsta(:),jlen(:),jvar(:)
     real(kind=8), allocatable :: c(:),g(:),jval(:),scaledc(:),scaledg(:), &
          scaledjval(:),xf(:),xg(:),xc(:),xj(:)
  end type memev_backup_type
  
  public :: memevini,memevend,memevscale,memevalf,memevalg,memevalc,memevalj, &
       memevinib,memevendb,memevbackup,memevrestore,memeviniconstr,memevendconstr, &
       memevscaleconstr
  
contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine memevini(evalf,evalg,evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)

    implicit none
    
    ! PROCEDURES
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    integer, intent(in) :: jnnzmax,m,n,p
    type(memev_type), intent(inout) :: memev
    type(c_ptr), optional, target, intent(in) :: pdata

    ! LOCAL SCALARS
    integer :: allocerr
    
    allocate(memev%jind(m+p),memev%jsorted(m+p),memev%jsta(m+p),memev%jlen(m+p), &
         memev%jvar(jnnzmax),memev%c(m+p),memev%g(n),memev%jval(jnnzmax), &
         memev%scaledc(m+p),memev%scaledg(n),memev%scaledjval(jnnzmax),memev%xf(n), &
         memev%xg(n),memev%xc(n),memev%xj(n),memev%sc(m+p),stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memevini, allocation error.'
       stop
    end if

    memev%evalf => evalf
    memev%evalg => evalg
       
    memev%evalc => evalc
    memev%evalj => evalj
    memev%evalhl => evalhl

    memev%jnnzmax = jnnzmax
    memev%n = n
    memev%m = m
    memev%p = p

    memev%pdatapresent = .false.
    if ( present( pdata ) ) then
       memev%pdatapresent = .true.
       memev%pdata => pdata
    end if

    memev%scale = .false.

    memev%firstf = .true.
    memev%firstg = .true.
    memev%firstc = .true.
    memev%firstj = .true.

  end subroutine memevini

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevend(memev)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(memev_type), intent(inout) :: memev

    ! LOCAL SCALARS
    integer :: allocerr
    
    deallocate(memev%jind,memev%jsorted,memev%jsta,memev%jlen,memev%jvar, &
         memev%c,memev%g,memev%jval,memev%scaledc,memev%scaledg, &
         memev%scaledjval,memev%xf,memev%xg,memev%xc,memev%xj,memev%sc, &
         stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memevend, deallocation error.'
       stop
    end if
    
  end subroutine memevend

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevscale(memev,n,x,m,p,ierr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: a,b,j
    
    ! LOCAL ARRAYS
    logical :: workl(m+p)

    if ( m + p .eq. 0 ) then
       write(*,*) 'There is no scaling for unconstrained problems.'
       return
    end if
    
    memev%scale = .false.

    call memevalg(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

!!$    memev%sf = max( 1.0d-8, 1.0d0 / max( 1.0d0, maxval( abs( memev%g(1:n) ) ) ) )

    memev%sf = max( 1.0d-8, min( 100.0d0 / maxval( abs( memev%g(1:n) ) ), 1.0d0 ) ) 

    memev%scaledg(1:n) = memev%g(1:n) * memev%sf
    
    workl(1:m+p) = .true.
    call memevalj(n,x,m,p,workl,ierr,memev)
    if ( ierr .ne. 0 ) return

    do j = 1,m + p
       a = memev%jsta(j)
       b = a + memev%jlen(j) - 1
!!$       memev%sc(j) = max( 1.0d-8, 1.0d0 / max( 1.0d0, maxval( abs( memev%jval(a:b) ) ) ) )
       
       memev%sc(j) = max( 1.0d-8, min( 100.0d0 / maxval( abs( memev%jval(a:b) ) ), 1.0d0 ) )
       
       memev%scaledjval(a:b) = memev%jval(a:b) * memev%sc(j)
    end do

    memev%scale = .true.
       
    write(*,*) 'Scaling factors were computed.'
    write(*,*) 'sf = ',memev%sf
    write(*,*) 'min_{j=1,m+p} { sc_j } = ',minval( memev%sc(1:m+p) )
    
  end subroutine memevscale

  ! *****************************************************************
  ! *****************************************************************

  subroutine memeviniconstr(evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)

    implicit none
    
    ! PROCEDURES
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    integer, intent(in) :: jnnzmax,m,n,p
    type(memev_type), intent(inout) :: memev
    type(c_ptr), optional, target, intent(in) :: pdata

    ! LOCAL SCALARS
    integer :: allocerr
    
    allocate(memev%jind(m+p),memev%jsorted(m+p),memev%jsta(m+p),memev%jlen(m+p), &
         memev%jvar(jnnzmax),memev%c(m+p),memev%jval(jnnzmax),memev%scaledc(m+p), &
         memev%scaledjval(jnnzmax),memev%xc(n),memev%xj(n),memev%sc(m+p),stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memeviniconstr, allocation error.'
       stop
    end if

    memev%evalc => evalc
    memev%evalj => evalj
    memev%evalhl => evalhl

    memev%jnnzmax = jnnzmax
    memev%n = n
    memev%m = m
    memev%p = p

    memev%pdatapresent = .false.
    if ( present( pdata ) ) then
       memev%pdatapresent = .true.
       memev%pdata => pdata
    end if

    memev%scale = .false.

    memev%firstc = .true.
    memev%firstj = .true.

  end subroutine memeviniconstr

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevendconstr(memev)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(memev_type), intent(inout) :: memev

    ! LOCAL SCALARS
    integer :: allocerr
    
    deallocate(memev%jind,memev%jsorted,memev%jsta,memev%jlen,memev%jvar, &
         memev%c,memev%jval,memev%scaledc,memev%scaledjval,memev%xc,memev%xj, &
         memev%sc,stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memevendconstr, deallocation error.'
       stop
    end if
    
  end subroutine memevendconstr

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevscaleconstr(memev,n,x,m,p,ierr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: a,b,j
    
    ! LOCAL ARRAYS
    logical :: workl(m+p)

    if ( m + p .eq. 0 ) then
       write(*,*) 'There is no scaling for unconstrained problems.'
       return
    end if
    
    memev%scale = .false.

    workl(1:m+p) = .true.
    call memevalj(n,x,m,p,workl,ierr,memev)
    if ( ierr .ne. 0 ) return

    do j = 1,m + p
       a = memev%jsta(j)
       b = a + memev%jlen(j) - 1
!!$       memev%sc(j) = max( 1.0d-8, 1.0d0 / max( 1.0d0, maxval( abs( memev%jval(a:b) ) ) ) )
       
       memev%sc(j) = max( 1.0d-8, min( 100.0d0 / maxval( abs( memev%jval(a:b) ) ), 1.0d0 ) )

       memev%scaledjval(a:b) = memev%jval(a:b) * memev%sc(j)
    end do

    memev%scale = .true.
       
    write(*,*) 'Scaling factors for the constraints were computed.'
    write(*,*) 'min_{j=1,m+p} { sc_j } = ',minval( memev%sc(1:m+p) )
    
  end subroutine memevscaleconstr

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevalf(n,x,ierr,memev)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    logical :: xdiff

    xdiff = .true.
    if ( .not. memev%firstf ) then
       if ( all( x(1:n) .eq. memev%xf(1:n) ) ) then
          xdiff = .false.
       end if
    end if
    
    if ( xdiff ) then
       memev%firstf = .false.
       memev%xf(1:n) = x(1:n)

       if ( memev%pdatapresent ) then
          call memev%evalf(n,x,memev%f,ierr,memev%pdata)
       else
          call memev%evalf(n,x,memev%f,ierr)
       end if

       if ( memev%scale ) then
          memev%scaledf = memev%f * memev%sf
       else
          memev%scaledf = memev%f
       end if
    end if

  end subroutine memevalf
 
  ! *****************************************************************
  ! *****************************************************************

  subroutine memevalg(n,x,ierr,memev)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    logical :: xdiff

    xdiff = .true.
    if ( .not. memev%firstg ) then
       if ( all( x(1:n) .eq. memev%xg(1:n) ) ) then
          xdiff = .false.
       end if
    end if
    
    if ( xdiff ) then
       memev%firstg = .false.
       memev%xg(1:n) = x(1:n)
       
       if ( memev%pdatapresent ) then
          call memev%evalg(n,x,memev%g,ierr,memev%pdata)
       else
          call memev%evalg(n,x,memev%g,ierr)
       end if

       if ( memev%scale ) then
          memev%scaledg(1:n) = memev%g(1:n) * memev%sf
       else
          memev%scaledg(1:n) = memev%g(1:n)
       end if
    end if

  end subroutine memevalg
 
  ! *****************************************************************
  ! *****************************************************************

  subroutine memevalc(n,x,ierr,memev)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    logical :: xdiff
    integer :: m,p
    
    xdiff = .true.
    if ( .not. memev%firstc ) then
       if ( all( x(1:n) .eq. memev%xc(1:n) ) ) then
          xdiff = .false.
       end if
    end if
    
    if ( xdiff ) then
       memev%firstc = .false.
       memev%xc(1:n) = x(1:n)
       
       m = memev%m
       p = memev%p
       
       if ( memev%pdatapresent ) then
          call memev%evalc(n,x,m,p,memev%c,ierr,memev%pdata)
       else
          call memev%evalc(n,x,m,p,memev%c,ierr)
       end if

       if ( memev%scale ) then
          memev%scaledc(1:m+p) = memev%c(1:m+p) * memev%sc(1:m+p)
       else
          memev%scaledc(1:m+p) = memev%c(1:m+p)
       end if
    end if

  end subroutine memevalc
 
  ! *****************************************************************
  ! *****************************************************************

  subroutine memevalj(n,x,m,p,jind,ierr,memev)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: m,n,p
    integer, intent(inout) :: ierr
    type(memev_type), intent(inout) :: memev
    
    ! ARRAY ARGUMENTS
    logical, intent(in) :: jind(m+p)
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    logical :: xdiff
    integer :: a,b,j
    
    xdiff = .true.
    if ( .not. memev%firstj ) then
       if ( all( x(1:n) .eq. memev%xj(1:n) ) .and. all( .not. jind(1:m+p) .or. memev%jind(1:m+p) ) ) then
          xdiff = .false.
       end if
    end if
    
    if ( xdiff ) then
       memev%firstj = .false.
       memev%xj(1:n) = x(1:n)
       memev%jind(1:m+p) = jind(1:m+p)

       if ( memev%pdatapresent ) then
          call memev%evalj(n,x,m,p,memev%jind,memev%jsorted,memev%jsta, &
               memev%jlen,memev%jnnzmax,memev%jvar,memev%jval,ierr,memev%pdata)
       else
          call memev%evalj(n,x,m,p,memev%jind,memev%jsorted,memev%jsta, &
               memev%jlen,memev%jnnzmax,memev%jvar,memev%jval,ierr)
       end if

       do j = 1,m + p
          if ( memev%jind(j) ) then
             a = memev%jsta(j)
             b = a + memev%jlen(j) - 1

             if ( memev%scale ) then
                memev%scaledjval(a:b) = memev%jval(a:b) * memev%sc(j)
             else
                memev%scaledjval(a:b) = memev%jval(a:b)
             end if
          end if
       end do
    end if

  end subroutine memevalj

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevinib(memev,memevb)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(memev_type), intent(in) :: memev
    type(memev_backup_type), intent(out) :: memevb

    ! LOCAL SCALARS
    integer :: allocerr,jnnzmax,m,n,p

    n = memev%n
    m = memev%m
    p = memev%p
    jnnzmax = memev%jnnzmax
    
    allocate(memevb%jind(m+p),memevb%jsorted(m+p),memevb%jsta(m+p),memevb%jlen(m+p), &
         memevb%jvar(jnnzmax),memevb%c(m+p),memevb%g(n),memevb%jval(jnnzmax), &
         memevb%scaledc(m+p),memevb%scaledg(n),memevb%scaledjval(jnnzmax),memevb%xf(n), &
         memevb%xg(n),memevb%xc(n),memevb%xj(n),stat=allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memevinib, allocation error.'
       stop
    end if
    
  end subroutine memevinib

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevendb(memevb)

    implicit none
    
    ! SCALAR ARGUMENTS
    type(memev_backup_type), intent(inout) :: memevb

    ! LOCAL SCALARS
    integer :: allocerr

    deallocate(memevb%jind,memevb%jsorted,memevb%jsta,memevb%jlen, &
         memevb%jvar,memevb%c,memevb%g,memevb%jval,memevb%scaledc, &
         memevb%scaledg,memevb%scaledjval,memevb%xf,memevb%xg, &
         memevb%xc,memevb%xj,stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In memevendb, deallocation error.'
       stop
    end if
    
  end subroutine memevendb

  ! *****************************************************************
  ! *****************************************************************

  subroutine memevbackup(memev,memevb)

    implicit none

    ! SCALARS ARGUMENTS
    type(memev_type), intent(in) :: memev
    type(memev_backup_type), intent(inout) :: memevb

    ! LOCAL SCALARS
    integer :: a,b,j,m,n,p

    n = memev%n
    m = memev%m
    p = memev%p
    
    memevb%firstf = memev%firstf
    if ( .not. memevb%firstf ) then
       memevb%xf(1:n) = memev%xf(1:n)
       memevb%f = memev%f
       memevb%scaledf = memev%scaledf
    end if

    memevb%firstg = memev%firstg
    if ( .not. memevb%firstg ) then
       memevb%xg(1:n) = memev%xg(1:n)
       memevb%g(1:n) = memev%g(1:n)
       memevb%scaledg(1:n) = memev%scaledg(1:n)
    end if
    
    memevb%firstc = memev%firstc
    if ( .not. memevb%firstc ) then
       memevb%xc(1:n) = memev%xc(1:n)
       memevb%c(1:m+p) = memev%c(1:m+p)
       memevb%scaledc(1:m+p) = memev%scaledc(1:m+p)
    end if
    
    memevb%firstj = memev%firstj
    if ( .not. memevb%firstj ) then
       memevb%xj(1:n) = memev%xj(1:n)
       memevb%jind(1:m+p) = memev%jind(1:m+p)
       do j = 1,m + p
          if ( memevb%jind(j) ) then
             memevb%jsorted(j) = memev%jsorted(j)
             memevb%jsta(j) = memev%jsta(j)
             memevb%jlen(j) = memev%jlen(j)

             a = memevb%jsta(j)
             b = a + memevb%jlen(j) - 1
             memevb%jvar(a:b) = memev%jvar(a:b)
             memevb%jval(a:b) = memev%jval(a:b)
             memevb%scaledjval(a:b) = memev%scaledjval(a:b)
          end if
       end do
    end if
    
  end subroutine memevbackup
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine memevrestore(memev,memevb)

    implicit none

    ! SCALARS ARGUMENTS
    type(memev_type), intent(inout) :: memev
    type(memev_backup_type), intent(in) :: memevb

    ! LOCAL SCALARS
    integer :: a,b,j,m,n,p

    n = memev%n
    m = memev%m
    p = memev%p
    
    memev%firstf = memevb%firstf
    if ( .not. memev%firstf ) then
       memev%xf(1:n) = memevb%xf(1:n)
       memev%f = memevb%f
       memev%scaledf = memevb%scaledf
    end if

    memev%firstg = memevb%firstg
    if ( .not. memev%firstg ) then
       memev%xg(1:n) = memevb%xg(1:n)
       memev%g(1:n) = memevb%g(1:n)
       memev%scaledg(1:n) = memevb%scaledg(1:n)
    end if
    
    memev%firstc = memevb%firstc
    if ( .not. memev%firstc ) then
       memev%xc(1:n) = memevb%xc(1:n)
       memev%c(1:m+p) = memevb%c(1:m+p)
       memev%scaledc(1:m+p) = memevb%scaledc(1:m+p)
    end if
    
    memev%firstj = memevb%firstj
    if ( .not. memev%firstj ) then
       memev%xj(1:n) = memevb%xj(1:n)
       memev%jind(1:m+p) = memevb%jind(1:m+p)
       do j = 1,m + p
          if ( memev%jind(j) ) then
             memev%jsorted(j) = memevb%jsorted(j)
             memev%jsta(j) = memevb%jsta(j)
             memev%jlen(j) = memevb%jlen(j)

             a = memev%jsta(j)
             b = a + memev%jlen(j) - 1
             memev%jvar(a:b) = memevb%jvar(a:b)
             memev%jval(a:b) = memevb%jval(a:b)
             memev%scaledjval(a:b) = memevb%scaledjval(a:b)
          end if
       end do
    end if
    
  end subroutine memevrestore
  
end module memev

