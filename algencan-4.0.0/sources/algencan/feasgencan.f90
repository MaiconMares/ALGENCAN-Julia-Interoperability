module bmfeasgencan

  use lss, only: lss_type,lssini,lssend
  use memev, only: memev_type,memeviniconstr,memevendconstr,memevscaleconstr,memevalc,memevalj
  use bmgencan, only: gencanb,genuncb
  use iso_c_binding, only: c_ptr,c_loc,c_f_pointer
  
  implicit none

  interface
     subroutine fsub(n,x,f,ierr,pdata) bind(C)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: f 
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine fsub

     subroutine gsub(n,x,g,ierr,pdata) bind(C)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: g(n)
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine gsub

     subroutine csub(n,x,m,p,c,ierr,pdata) bind(C)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: m,n,p
       real(kind=8), intent(in) :: x(n)
       real(kind=8), intent(out) :: c(m+p) 
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine csub

     subroutine jsub(n,x,m,p,ind,jsorted,jsta,jlen,lim,jvar,jval,ierr,pdata) bind(C)
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

     subroutine hlsub(n,x,m,p,lambda,lim,inclf,hlnnz,hlrow,hlcol,hlval,ierr,pdata) bind(C)
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

  type, private :: feas_type
     logical :: avoidhsqf = .false.
     integer :: nmin = 200
     real(kind=8) :: jdenmax = 0.2d0, hdenmax = 0.1d0
     
     real(kind=8) :: epsfeas
     logical, dimension(:), pointer :: workl
     real(kind=8), dimension(:), pointer :: workr
     type(memev_type), pointer :: memev
  end type feas_type
  
  private
  
  ! SCALAR PARAMETERS
  integer, parameter :: lrgstint = huge( 1 )
  real(kind=8), parameter :: lrgstreal = huge( 1.0d0 )
  
  public :: feasgencan,feasgenunc,feasgencanb,feasgenuncb

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine feasgenunc(evalc,evalj,evalhl,jnnzmax,hlnnzmax,n,x,m,p, &
       epsfeas,maxit,scale,extallowed,csupn,sqf,nsqf,nsqfpsupn,iter, &
       ierr,istop,pdata) bind(C, name="feasgenunc")

    implicit none
    
    ! PROCEDURES
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,scale
    integer, intent(in) :: hlnnzmax,jnnzmax,m,maxit,n,p
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: epsfeas
    real(kind=8), intent(out) :: csupn,nsqfpsupn,sqf
    type(c_ptr), optional, target, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: nsqf(n)

    ! LOCAL SCALARS
    type(memev_type) :: memev
    type(lss_type) :: lss

    ierr = 0

    ! This wrong, since more space is needed
    call lssini(hlnnzmax,lss)

    call memeviniconstr(evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)    

    if ( scale ) then
       call memevscaleconstr(memev,n,x,m,p,ierr)
       if ( ierr .ne. 0 ) return
    end if

    call feasgenuncb(lss,memev,n,x,m,p,epsfeas,maxit, &
         extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)
    
    if ( ierr .ne. 0 ) return

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    csupn = max( maxval( abs( memev%c(1:m) ) ), maxval( max( 0.0d0, memev%c(m+1:m+p) ) ) )
    
    call memevendconstr(memev)
    
    call lssend(lss)
    
  end subroutine feasgenunc
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine feasgenuncb(lss,memev,n,x,m,p,epsfeas,maxit, &
       extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed
    integer, intent(in) :: m,maxit,n,p
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: epsfeas
    real(kind=8), intent(out) :: nsqfpsupn,sqf
    type(memev_type), intent(inout), target :: memev
    type(lss_type), intent(inout) :: lss

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: nsqf(n)

    ! LOCAL SCALARS
    logical :: hfixstr
    real(kind=8) :: epsopt,ftarget
    type(feas_type), target :: feas

    ! LOCAL ARRAYS
    logical, target :: workl(m+p)
    real(kind=8), target :: workr(m+p)

    feas%epsfeas = epsfeas
    feas%workr => workr(1:m+p)
    feas%workl => workl(1:m+p)
    feas%memev => memev
    
    hfixstr = .false.
    epsopt = 0.0d0
    ftarget = 0.0d0

    call genuncb(evalsqf,evalnsqf,evalhsqf,hfixstr,n,x,sqf,nsqf,nsqfpsupn, &
         lss,ftarget,epsopt,maxit,extallowed,iter,ierr,istop,evalstp,c_loc(feas))
    
  end subroutine feasgenuncb
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine feasgencan(evalc,evalj,evalhl,jnnzmax,hlnnzmax,n,x,lind, &
       lbnd,uind,ubnd,m,p,epsfeas,maxit,scale,extallowed,csupn,sqf,nsqf, &
       nsqfpsupn,iter,ierr,istop,pdata) bind(C, name="feasgencan")

    implicit none
    
    ! PROCEDURES
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,scale
    integer, intent(in) :: hlnnzmax,jnnzmax,m,maxit,n,p
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: epsfeas
    real(kind=8), intent(out) :: csupn,nsqfpsupn,sqf
    type(c_ptr), optional, target, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: nsqf(n)

    ! LOCAL SCALARS
    type(memev_type) :: memev
    type(lss_type) :: lss

    ierr = 0

    ! This wrong, since more space is needed
    call lssini(hlnnzmax,lss)

    call memeviniconstr(evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)    

    if ( scale ) then
       call memevscaleconstr(memev,n,x,m,p,ierr)
       if ( ierr .ne. 0 ) return
    end if

    call feasgencanb(lss,memev,n,x,lind,lbnd,uind,ubnd,m,p,epsfeas, &
       maxit,extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)
    
    if ( ierr .ne. 0 ) return

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    csupn = max( maxval( abs( memev%c(1:m) ) ), maxval( max( 0.0d0, memev%c(m+1:m+p) ) ) )
    
    call memevendconstr(memev)
    
    call lssend(lss)
    
  end subroutine feasgencan
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine feasgencanb(lss,memev,n,x,lind,lbnd,uind,ubnd,m,p,epsfeas, &
       maxit,extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)

    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed
    integer, intent(in) :: m,maxit,n,p
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: epsfeas
    real(kind=8), intent(out) :: nsqfpsupn,sqf
    type(memev_type), intent(inout), target :: memev
    type(lss_type), intent(inout) :: lss

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: nsqf(n)

    ! LOCAL SCALARS
    logical :: hfixstr
    real(kind=8) :: epsopt,ftarget
    type(feas_type), target :: feas

    ! LOCAL ARRAYS
    logical, target :: workl(m+p)
    real(kind=8), target :: workr(m+p)

    feas%epsfeas = epsfeas
    feas%workr => workr(1:m+p)
    feas%workl => workl(1:m+p)
    feas%memev => memev
    
    hfixstr = .false.
    epsopt = 0.0d0
    ftarget = 0.0d0

    call gencanb(evalsqf,evalnsqf,evalhsqf,hfixstr,n,x,lind,lbnd,uind,ubnd, &
         sqf,nsqf,nsqfpsupn,lss,ftarget,epsopt,maxit,extallowed,iter,ierr, &
         istop,evalstp,c_loc(feas))
    
  end subroutine feasgencanb
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalstp(n,x,gsupn,inhdefstp,stp,ierr,feasptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    real(kind=8), intent(in) :: gsupn
    logical, intent(out) :: inhdefstp,stp 
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: feasptr
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: m,p
    real(kind=8) :: csupn
    type(feas_type), pointer :: feas
    type(memev_type), pointer :: memev
    
    call c_f_pointer(feasptr,feas)
    memev => feas%memev
    
    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    csupn = max( maxval( abs( memev%c(1:m) ) ), maxval( max( 0.0d0, memev%c(m+1:m+p) ) ) )

    inhdefstp = .false.
    
    if ( csupn .le. feas%epsfeas .or. ( csupn .gt. sqrt( feas%epsfeas ) .and. &
         gsupn .le. feas%epsfeas * sqrt( feas%epsfeas ) ) ) then
       stp = .true.
    else
       stp = .false.
    end if
    
    open(20,file='solver-interrupted-tabline.txt')
    write(20,9000) n,m+p,csupn
    close(20)
  
9000 format(2(1X,I6),1X,1P,D7.1)

  end subroutine evalstp

  ! *****************************************************************
  ! *****************************************************************

  subroutine evalsqf(n,x,sqf,ierr,feasptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    real(kind=8), intent(out) :: sqf
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: feasptr
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: m,p
    type(feas_type), pointer :: feas
    type(memev_type), pointer :: memev

    if ( .not. present( feasptr ) ) then
       write(*,*) 'There is something wrong in evalsqf.'
       stop
    end if
    
    call c_f_pointer(feasptr,feas)
    memev => feas%memev
    
    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    sqf = 0.5d0 * ( sum( memev%c(1:m) ** 2 ) + sum( max( 0.0d0, memev%c(m+1:m+p) ) ** 2 ) )

  end subroutine evalsqf
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalnsqf(n,x,nsqf,ierr,feasptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: feasptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: nsqf(n)

    ! LOCAL SCALARS
    integer :: a,b,j,m,p
    type(feas_type), pointer :: feas
    type(memev_type), pointer :: memev

    if ( .not. present( feasptr ) ) then
       write(*,*) 'There is something wrong in evalnsqf.'
       stop
    end if
    
    call c_f_pointer(feasptr,feas)
    memev => feas%memev

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    feas%workr(1:m) = memev%c(1:m)
    feas%workr(m+1:m+p) = max( 0.0d0, memev%c(m+1:m+p) )

    feas%workl(1:m+p) = feas%workr(1:m+p) .ne. 0.0d0

    nsqf(1:n) = 0.0d0

    if ( any( feas%workl(1:m+p) ) ) then
       call memevalj(n,x,m,p,feas%workl,ierr,memev)
       if ( ierr .ne. 0 ) return

       do j = 1,m + p
          if ( feas%workr(j) .ne. 0.0d0 ) then
             a = memev%jsta(j)
             b = a + memev%jlen(j) - 1
             nsqf(memev%jvar(a:b)) = nsqf(memev%jvar(a:b)) + feas%workr(j) * memev%jval(a:b)
          end if
       end do
    end if

  end subroutine evalnsqf
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhsqf(n,x,lim,hsqfnnz,hsqfrow,hsqfcol,hsqfval,ierr,feasptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,n
    integer, intent(inout) :: ierr
    integer, intent(out) :: hsqfnnz
    type(c_ptr), optional, intent(in) :: feasptr
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    integer, intent(out) :: hsqfcol(lim),hsqfrow(lim)
    real(kind=8), intent(out) :: hsqfval(lim)

    ! LOCAL SCALARS
    integer :: a,b,dim,hcnnz,j,jdensit,jlenmax,k,l,m,nelem,p
    type(feas_type), pointer :: feas
    type(memev_type), pointer :: memev

    if ( .not. present( feasptr ) ) then
       write(*,*) 'There is something wrong in evalhsqf.'
       stop
    end if
    
    call c_f_pointer(feasptr,feas)
    memev => feas%memev

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    feas%workr(1:m) = memev%c(1:m)
    feas%workr(m+1:m+p) = max( 0.0d0, memev%c(m+1:m+p) )

    if ( memev%pdatapresent ) then
       call memev%evalhl(n,x,m,p,feas%workr,lim,.false.,hcnnz,hsqfrow,hsqfcol,hsqfval,ierr, &
            memev%pdata)
    else
       call memev%evalhl(n,x,m,p,feas%workr,lim,.false.,hcnnz,hsqfrow,hsqfcol,hsqfval,ierr)
    end if

    hsqfnnz = hcnnz
    
    if ( ierr .ne. 0 ) return
          
    feas%workl(1:m+p) = feas%workr(1:m+p) .ne. 0.0d0

    if ( any( feas%workl(1:m+p) ) ) then
       
       call memevalj(n,x,m,p,feas%workl,ierr,memev)
       if ( ierr .ne. 0 ) return

       jlenmax = maxval( memev%jlen(1:m+p), feas%workl(1:m+p) )
       jdensit = sum( ( memev%jlen(1:m+p) / n ) * ( ( memev%jlen(1:m+p) + 1 ) / ( n + 1 ) ), feas%workl(1:m+p) )
       
       if ( .not. feas%avoidhsqf .and. ( n .le. feas%nmin .or. &
            ( jlenmax .le. feas%jdenmax * n .and. jdensit .le. feas%hdenmax ) ) ) then

          hsqfnnz = hcnnz
          
          outer: do j = 1,m + p
             if ( feas%workl(j) ) then
                a = memev%jsta(j)
                b = a + memev%jlen(j) - 1
                do k = a,b
                   if ( memev%jval(k) .ne. 0.0d0 ) then
                      if ( .not. memev%jsorted(j) ) then
                         do l = a,b
                            if ( memev%jvar(l) .ge. memev%jvar(k) .and. memev%jval(l) .ne. 0.0d0 ) then
!!$                               if (  hsqfnnz + 1 .gt. lim .or. .not. ( n .le. feas%nmin .or. &
!!$                                    ( 2.0d0 * ( hsqfnnz + 1 - hcnnz ) ) / n / ( n + 1 ) .le. feas%hdenmax ) ) then
                               if (  hsqfnnz + 1 .gt. lim ) then
                                  feas%avoidhsqf = .true.
                                  exit outer
                               end if
                               
                               hsqfnnz = hsqfnnz + 1
                               hsqfval(hsqfnnz) = memev%jval(k) * memev%jval(l)
                               hsqfcol(hsqfnnz) = memev%jvar(l)
                               hsqfrow(hsqfnnz) = memev%jvar(k)
                            end if
                         end do
                      else
                         do l = k,b
                            if ( memev%jval(l) .ne. 0.0d0 ) then
!!$                                if (  hsqfnnz + 1 .gt. lim .or. .not. ( n .le. feas%nmin .or. &
!!$                                    ( 2.0d0 * ( hsqfnnz + 1 - hcnnz ) ) / n / ( n + 1 ) .le. feas%hdenmax ) ) then
                                if (  hsqfnnz + 1 .gt. lim ) then
                                  feas%avoidhsqf = .true.
                                  exit outer
                               end if
                               
                               hsqfnnz = hsqfnnz + 1
                               hsqfval(hsqfnnz) = memev%jval(k) * memev%jval(l)
                               hsqfcol(hsqfnnz) = memev%jvar(l)
                               hsqfrow(hsqfnnz) = memev%jvar(k)
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end do outer

          if ( .not. feas%avoidhsqf ) then
             call collapse(n,hsqfnnz,hsqfrow,hsqfcol,hsqfval)
             call removenull(hsqfnnz,hsqfrow,hsqfcol,hsqfval)
          end if
       end if

       if ( feas%avoidhsqf .or. .not. ( n .le. feas%nmin .or. &
            ( jlenmax .le. feas%jdenmax * n .and. jdensit .le. feas%hdenmax ) ) ) then
          
          dim = n
    
          hsqfnnz = hcnnz
          
          do j = 1,m + p
             if ( feas%workl(j) ) then
                nelem = memev%jlen(j)
             
                if ( hsqfnnz + nelem + 1 .gt. lim ) then
                   ierr = -98
                   return
                end if
             
                dim = dim + 1

                a = memev%jsta(j)
                b = a + nelem - 1

                hsqfrow(hsqfnnz+1:hsqfnnz+nelem) = memev%jvar(a:b)
                hsqfcol(hsqfnnz+1:hsqfnnz+nelem) = dim
                hsqfval(hsqfnnz+1:hsqfnnz+nelem) = memev%jval(a:b)

                hsqfrow(hsqfnnz+nelem+1) = dim
                hsqfcol(hsqfnnz+nelem+1) = dim
                hsqfval(hsqfnnz+nelem+1) = - 1.0d0
                
                hsqfnnz = hsqfnnz + nelem + 1
             end if
          end do

       end if

    end if

  end subroutine evalhsqf

  ! *****************************************************************
  ! *****************************************************************

  subroutine collapse(n,hsqfnnz,hsqfrow,hsqfcol,hsqfval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: hsqfnnz,n

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hsqfcol(hsqfnnz),hsqfrow(hsqfnnz)
    real(kind=8), intent(inout) :: hsqfval(hsqfnnz)

    ! LOCAL SCALARS
    integer :: col,i,k,next,row
    
    ! LOCAL ARRAYS
    integer :: hsqfnext(hsqfnnz),hsqfrowsta(n)
    real(kind=8) :: r(n)
    
    hsqfrowsta(1:n) = 0
       
    do k = 1,hsqfnnz
       row = hsqfrow(k)
       hsqfnext(k) = hsqfrowsta(row)
       hsqfrowsta(row) = k
    end do

    r(1:n) = 0.0d0
    
    do i = 1,n
       next = hsqfrowsta(i)
       do while ( next .ne. 0 )
          col = hsqfcol(next)
          r(col) = r(col) + hsqfval(next)
          next = hsqfnext(next)
       end do

       next = hsqfrowsta(i)
       do while ( next .ne. 0 )
          col = hsqfcol(next)
          hsqfval(next) = r(col)
          r(col) = 0.0d0
          next = hsqfnext(next)
       end do
    end do
       
  end subroutine collapse

  ! *****************************************************************
  ! *****************************************************************

  subroutine removenull(hsqfnnz,hsqfrow,hsqfcol,hsqfval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: hsqfnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: hsqfcol(hsqfnnz),hsqfrow(hsqfnnz)
    real(kind=8), intent(inout) :: hsqfval(hsqfnnz)

    ! LOCAL SCALARS
    integer :: k,tmp
    
    ! LOCAL ARRAYS
    logical :: nonnull(hsqfnnz)
    integer :: nnind(hsqfnnz)

    nonnull(1:hsqfnnz) = hsqfval(1:hsqfnnz) .ne. 0.0d0
    
    tmp = count( nonnull(1:hsqfnnz) )
    nnind(1:tmp) = pack( (/ (k, k=1,hsqfnnz) /), nonnull(1:hsqfnnz) )
       
    hsqfrow(1:tmp) = hsqfrow(nnind(1:tmp))
    hsqfcol(1:tmp) = hsqfcol(nnind(1:tmp))
    hsqfval(1:tmp) = hsqfval(nnind(1:tmp))
    hsqfnnz = tmp

  end subroutine removenull

end module bmfeasgencan
