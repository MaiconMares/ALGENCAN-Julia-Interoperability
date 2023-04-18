module bmalgencan

  use lss, only: lss_type,lssini,lssend
  use memev, only: memev_type,memevini,memevend,memevscale,memevalf,memevalg,memevalc,memevalj, &
       memev_backup_type,memevinib,memevendb,memevbackup,memevrestore
  use bmnewtkkt, only: newtkktb
  use bmgencan, only: gencanb,project
  use bmfeasgencan, only: feasgencanb,feasgenuncb
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

  type, private :: aldata_type
     logical :: avoidhal = .false.
     integer :: nmin = 200
     real(kind=8) :: jdenmax = 0.2d0, hdenmax = 0.1d0
     
     real(kind=8), pointer :: rho
     logical, dimension(:), pointer :: workl
     real(kind=8), dimension(:), pointer :: lambda,workr
     type(memev_type), pointer :: memev
  end type aldata_type
  
  private
  
  ! SCALAR PARAMETERS
  integer, parameter :: lrgstint = huge( 1 )
  real(kind=8), parameter :: lrgstreal = huge( 1.0d0 ), macheps = epsilon( 1.0d0 ), &
       macheps12 = sqrt( macheps ), machepsinv = 1.0d0 / macheps,  &
       macheps12inv = 1.0d0 / macheps12, lammin = - machepsinv, lammax = machepsinv, &
       rhofrac = 0.5d0, rhomult = 10.0d0, rhoinimin = macheps12, rhoinimax = macheps12inv
  
  public :: algencan

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine algencan(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
       n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt,maxoutit, &
       scale,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn, &
       bdsvio,outiter,totiter,nwcalls,nwtotit,ierr,istop,pdata) bind(C, name="algencan")

    implicit none
    
    ! PROCEDURES
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    logical, intent(in) :: corrin,extallowed,rhoauto,scale
    integer, intent(in) :: hlnnzmax,jnnzmax,m,maxoutit,n,p
    integer, intent(out) :: ierr,istop,nwcalls,nwtotit,outiter,totiter
    real(kind=8), intent(in) :: epsfeas,epscompl,epsopt,rhoini
    real(kind=8), intent(out) :: bdsvio,csupn,f,nlpsupn,ssupn
    type(c_ptr), optional, target, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: lambda(m+p),x(n)

    ! istop = 0 means an AL iterate satisfied KKT conditions

    ! istop = 1 means the acceleration process found a KKT point when
    ! an AL iterate satisfying KKT conditions with half precision has
    ! already been found

    ! istop = 2 means too many consecutive fails of the inner solver,
    ! but a KKT or a feasible point has been found
    
    ! istop = 3 means the maximum number of outer iterations was
    ! reached, but a KKT or a feasible point has been found

    ! istop = 4 means the penalty parameter is larger than 1e+20, but
    ! a KKT or a feasible point has been found

    ! Note: in cases 1, 2, and 3, if a KKT point was found then it was
    ! found by the acceleration process when the AL iterates were far
    ! from satisfying KKT conditions (less than half of the required
    ! precision). In this case, since stopping returning this KKT
    ! point was considered premature, the minimization process
    ! continued. At the end, it would have been more efficient to stop
    ! returning it, but how could we know in advance?

    ! If the AL approach stops without returning a feasible point, all
    ! it did is discarded; and the squared infeasibility is minimized
    ! starting from the given initial guess.

    ! istop = 5 means that the initial guess is in fact feasible (and
    ! that feasibility was lost by the AL approach)

    ! istop = 6 means the squared infeasibility is in fact
    ! minimized. Whether the final point is feasible or not depends on
    ! the reported feasibility, but istop = 6 in any case. The
    ! objective function at the final iterate of the
    ! squared-feasibility minimization process is evaluated and
    ! reported.
    
    ! LOCAL SCALARS
    logical :: kktfound,feasfound
    integer :: iter,maxit,nbds
    real(kind=8) :: bdsvioini,bestbdsvio,bestcsupn,bestf,bestssupn,bestnlpsupn, &
         csupnini,fini,nsqfpsupn,sqf
    type(memev_type) :: memev
    type(memev_backup_type) :: memevb
    type(lss_type) :: lss

    ! LOCAL ARRAYS
    real(kind=8) :: bestx(n),bestlambda(m+p),nsqf(n),xini(n)

    ! integer :: hlnnz
    ! real(kind=8) :: g(n)
    ! logical :: sorted(m+p), ind(m+p)
    ! integer :: jsta(m+p),jlen(m+p),jvar(jnnzmax)
    ! real(kind=8) :: jval(jnnzmax)
    ! integer :: hlrow(hlnnzmax),hlcol(hlnnzmax)
    ! real(kind=8) :: hlval(hlnnzmax)
    ! print *, "PARAMETERS VALUES"
    ! print *, "jnnzmax = ", jnnzmax
    ! print *, "hlnnzmax = ", hlnnzmax
    ! print *, "n = ", n
    ! print *, "x = ", x
    ! print *, "lind = ", lind
    ! print *, "lbnd = ", lbnd
    ! print *, "uind = ", uind
    ! print *, "ubnd = ", ubnd
    ! print *, "m = ", m
    ! print *, "p = ", p
    ! print *, "lambda = ", lambda
    ! print *, "epsfeas = ", epsfeas
    ! print *, "epscompl = ", epscompl
    ! print *, "epsopt = ", epsopt
    ! print *, "maxoutit = ", maxoutit
    ! print *, "scale = ", scale
    ! print *, "rhoauto = ", rhoauto
    ! print *, "rhoini = ", rhoini
    ! print *, "extallowed = ", extallowed
    ! print *, "corrin = ", corrin

    ! call evalf(n,x,f,ierr,pdata)
    ! print *, "f = ", f

    ! call evalg(n,x,g,ierr,pdata)
    ! print *, "g = ", g

    ! call evalc(n,x,m,p,g,ierr,pdata)
    ! print *, "c = ", g

    ! ind(1) = .true.
    ! call evalj(n,x,m,p,ind,sorted,jsta,jlen,jnnzmax,jvar,jval,ierr,pdata)
    ! print *, "sorted = ", sorted
    ! print *, "jsta = ", jsta
    ! print *, "jlen = ", jlen
    ! print *, "jvar = ", jvar
    ! print *, "jval = ", jval
    
    ! call evalhl(n,x,m,p,lambda,hlnnzmax,.true.,hlnnz,hlrow,hlcol,hlval,ierr,pdata)
    ! print *, "hlnnz = ", hlnnz
    ! print *, "hlrow = ", hlrow(1:hlnnz)
    ! print *, "hlcol = ", hlcol(1:hlnnz)
    ! print *, "hlval = ", hlval(1:hlnnz)
        
    ! stop
    
    ierr = 0

    call lssini(hlnnzmax,lss)
    
    call memevini(evalf,evalg,evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)    

    call project(n,lind,lbnd,uind,ubnd,x)
    
    xini(1:n) = x(1:n)
    bdsvioini = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

    if ( scale ) then
       call memevscale(memev,n,x,m,p,ierr)
       if ( ierr .ne. 0 ) return
    end if

    call memevinib(memev,memevb)
    
    call memevalf(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    fini = memev%f
    
    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    csupnini  = max( 0.0d0, max( maxval( abs( memev%c(1:m) ) ), maxval( memev%c(m+1:m+p) ) ) )
    
    call algencanb(lss,memev,memevb,n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas, &
         epscompl,epsopt,maxoutit,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn, &
         bdsvio,kktfound,feasfound,bestx,bestlambda,bestf,bestcsupn,bestssupn,bestnlpsupn, &
         bestbdsvio,outiter,totiter,nwcalls,nwtotit,ierr,istop)

    if ( ierr .ne. 0 ) return

    if ( kktfound .or. feasfound ) then
       x(1:n) = bestx(1:n)
       lambda(1:m+p) = bestlambda(1:m+p)

       f = bestf
       csupn = bestcsupn
       ssupn = bestssupn
       nlpsupn = bestnlpsupn
       bdsvio = bestbdsvio

    else if ( max( bdsvioini, csupnini ) .le. epsfeas ) then
       x(1:n) = xini(1:n)
       lambda(1:m+p) = 0.0d0
       
       f = fini
       bdsvio = bdsvioini
       csupn = csupnini
       
       istop = 5

    else
       x(1:n) = xini(1:n)

       maxit = 50000

       nbds = count( lind(1:n) ) + count( uind(1:n) )

       if ( nbds .eq. 0 ) then
          call feasgenuncb(lss,memev,n,x,m,p,epsfeas,maxit, &
               extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)

       else
          call feasgencanb(lss,memev,n,x,lind,lbnd,uind,ubnd,m,p,epsfeas, &
               maxit,extallowed,sqf,nsqf,nsqfpsupn,iter,ierr,istop)
       end if

       totiter = totiter + iter
       
       lambda(1:m+p) = 0.0d0

       call memevalf(n,x,ierr,memev)
       if ( ierr .ne. 0 ) return
          
       f = memev%f

       call memevalc(n,x,ierr,memev)
       if ( ierr .ne. 0 ) return

       csupn = max( 0.0d0, max( maxval( abs( memev%c(1:m) ) ), maxval( memev%c(m+1:m+p) ) ) )
       ssupn = 0.0d0
       nlpsupn = nsqfpsupn
       bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

       istop = 6
    end if
    
    call memevendb(memevb)
    
    call memevend(memev)
    
    call lssend(lss)
    
  end subroutine algencan
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine algencanb(lss,memev,memevb,n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas, &
       epscompl,epsopt,maxoutit,rhoauto,rhoini,extallowed,corrin,f,csupn,ssupn,nlpsupn, &
       bdsvio,kktfound,feasfound,bestx,bestlambda,bestf,bestcsupn,bestssupn,bestnlpsupn, &
       bestbdsvio,outiter,totiter,nwcalls,nwtotit,ierr,istop)

    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: corrin,extallowed,rhoauto
    logical, intent(out) :: kktfound,feasfound
    integer, intent(in) :: m,maxoutit,n,p
    integer, intent(inout) :: ierr
    integer, intent(out) :: istop,nwcalls,nwtotit,outiter,totiter
    real(kind=8), intent(in) :: epsfeas,epscompl,epsopt,rhoini
    real(kind=8), intent(out) :: bestf,bestbdsvio,bestcsupn,bestnlpsupn,bestssupn,f,bdsvio, &
         csupn,nlpsupn,ssupn
    type(memev_type), intent(inout), target :: memev
    type(memev_backup_type), intent(inout) :: memevb
    type(lss_type), intent(inout) :: lss

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout), target :: lambda(m+p)
    real(kind=8), intent(out), target :: bestlambda(m+p)
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: bestx(n)

    ! LOCAL SCALARS
    logical :: hfixstr,infeasfound,nwfeasfound,nwinfeasfound
    integer :: a,b,cnoprog,innistop,inniter,j,maxit,nwierr,nwistop,nwiter,nwmaxit,subfail
    real(kind=8) :: al,csupnprev,epsopk,ftarget,nfeaspsupn,nwbdsvio,nwcsupn,nwf,nwnlpsupn, &
         nwssupn,ssupnprev,sumc,nwfeasf,nwfeasbdsvio,nwfeascsupn,nwfeasssupn,nwfeasnlpsupn, &
         nwinfeasf,nwinfeasbdsvio,nwinfeascsupn,nwinfeasssupn,nwinfeasnlpsupn,feasf,feasbdsvio, &
         feascsupn,feasssupn,feasnlpsupn,kktf,kktbdsvio,kktcsupn,kktssupn,kktnlpsupn,infeasf, &
         infeasbdsvio,infeascsupn,infeasssupn,infeasnlpsupn
    real(kind=8), target :: rho
    type(aldata_type), target :: aldata
    
    ! LOCAL ARRAYS
    logical, target :: workl(m+p)
    real(kind=8), target :: workr(m+p)
    real(kind=8) :: feaslambda(m+p),feasx(n),infeaslambda(m+p),infeasx(n),kktlambda(m+p), &
         kktx(n),nfeas(n),nl(n),nwx(n),nwl(m+p),nwfeasx(n),nwfeasl(m+p),nwinfeasx(n),nwinfeasl(m+p)

    outiter = 0
    totiter = 0
    nwcalls = 0
    nwtotit = 0
    cnoprog = 0
    subfail = 0

    kktfound = .false.
    feasfound = .false.
    infeasfound = .false.
    
    call project(n,lind,lbnd,uind,ubnd,x)
    
    bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )
    
    call memevalf(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
          
    f = memev%f

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    lambda(m+1:m+p) = max( 0.0d0, lambda(m+1:m+p) )

    csupn = max( 0.0d0, max( maxval( abs( memev%c(1:m) ) ), maxval( memev%c(m+1:m+p) ) ) )
    ssupn = max( 0.0d0, maxval( min( - memev%c(m+1:m+p), lambda(m+1:m+p) ) ) )

    call memevalg(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    nl(1:n) = memev%scaledg(1:n)

    workl(1:m) = .true.
    workl(m+1:m+p) = lambda(m+1:m+p) .gt. 0.0d0

    if ( any( workl(1:m+p) ) ) then
       call memevalj(n,x,m,p,workl,ierr,memev)
       if ( ierr .ne. 0 ) return
    
       do j = 1,m + p
          if ( lambda(j) .ne. 0.0d0 ) then
             a = memev%jsta(j)
             b = a + memev%jlen(j) - 1
             nl(memev%jvar(a:b)) = nl(memev%jvar(a:b)) + lambda(j) * memev%scaledjval(a:b)
          end if
       end do
    end if

    nl(1:n) = x(1:n) - nl(1:n)

    call project(n,lind,lbnd,uind,ubnd,nl)

    nl(1:n) = nl(1:n) - x(1:n)
    
    nlpsupn = maxval( abs ( nl(1:n) ) )

    write(*,*) 'outiter = ',outiter,' csupn = ',csupn,' bdsvio = ',bdsvio,' ssupn = ',ssupn,' nlpsupn = ',nlpsupn
    write(*,*) 'f = ',f
    
    call updatesol(epsfeas,n,x,m,p,lambda,f,csupn,bdsvio,ssupn,nlpsupn, &
         kktfound,kktf,kktx,kktlambda,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
         feasfound,feasf,feasx,feaslambda,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
         infeasfound,infeasf,infeasx,infeaslambda,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
         bestf,bestx,bestlambda,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn, &
         epscompl,epsopt,outiter,totiter,nwcalls,nwtotit)

    write(*,*) 'kkt:    ',kktfound,kktf,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn
    write(*,*) 'feas:   ',feasfound,feasf,feascsupn,feasbdsvio,feasssupn,feasnlpsupn
    write(*,*) 'infeas: ',infeasfound,infeasf,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn
    write(*,*) 'best:   ',bestf,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn
    
    if ( max( bdsvio, csupn ) .le. epsfeas .and. ssupn .le. epscompl .and. nlpsupn .le. epsopt ) then
       write(*,*) 'Success (with the initial guess)!'
       istop = 0
       return
    end if

100 continue
  
    if ( outiter .ge. maxoutit ) then
       write(*,*) 'Maximum number of outer iterations reached!'
       istop = 3
       return
    end if

    write(*,*) 'Trying Newton-KKT.'
       
    nwcalls = nwcalls + 1
          
    call memevbackup(memev,memevb)
          
    nwx(1:n) = x(1:n)
    nwl(1:m+p) = lambda(1:m+p)

    if ( memev%scale ) then
       nwl(1:m+p) = memev%sf * nwl(1:m+p) / memev%sc(1:m+p)
    end if

    nwierr = 0
    nwmaxit = 10
    
    call newtkktb(lss,memev,n,nwx,lind,lbnd,uind,ubnd,m,p,nwl,epsfeas,epscompl, &
         epsopt,nwmaxit,corrin,nwcsupn,nwssupn,nwnlpsupn,nwbdsvio,nwfeasfound,nwfeasx, &
         nwfeasl,nwfeascsupn,nwfeasssupn,nwfeasnlpsupn,nwfeasbdsvio,nwinfeasfound,nwinfeasx, &
         nwinfeasl,nwinfeascsupn,nwinfeasssupn,nwinfeasnlpsupn,nwinfeasbdsvio,nwiter, &
         nwierr,nwistop)

    nwtotit = nwtotit + nwiter
          
    write(*,*)'Newton-KKT nwierr = ',nwierr,' nwistop = ',nwistop

    if ( nwfeasfound ) then
       if ( memev%scale ) then
          nwfeasl(1:m+p) = nwfeasl(1:m+p) * memev%sc(1:m+p) / memev%sf
       end if
          
       call memevalf(n,nwfeasx,ierr,memev)
       if ( ierr .ne. 0 ) return

       nwfeasf = memev%f

       call updatesol(epsfeas,n,nwfeasx,m,p,nwfeasl,nwfeasf,nwfeascsupn,nwfeasbdsvio,nwfeasssupn, &
            nwfeasnlpsupn,kktfound,kktf,kktx,kktlambda,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
            feasfound,feasf,feasx,feaslambda,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
            infeasfound,infeasf,infeasx,infeaslambda,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
            bestf,bestx,bestlambda,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn, &
            epscompl,epsopt,outiter,totiter,nwcalls,nwtotit)

       write(*,*) 'kkt:    ',kktfound,kktf,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn
       write(*,*) 'feas:   ',feasfound,feasf,feascsupn,feasbdsvio,feasssupn,feasnlpsupn
       write(*,*) 'infeas: ',infeasfound,infeasf,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn
       write(*,*) 'best:   ',bestf,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn

       if ( max( nwfeasbdsvio, nwfeascsupn ) .le. epsfeas .and. nwfeasssupn .le. epscompl .and. &
            nwfeasnlpsupn .le. epsopt .and. max( bdsvio, csupn ) .le. sqrt( epsfeas ) .and. &
            ssupn .le. sqrt( epscompl ) .and. nlpsupn .le. sqrt( epsopt ) ) then
          write(*,*) 'Success (after a Newton-KKT acceleration)!'
          istop = 1
          return
       end if

    else if ( nwinfeasfound ) then
       if ( memev%scale ) then
          nwinfeasl(1:m+p) = nwinfeasl(1:m+p) * memev%sc(1:m+p) / memev%sf
       end if
          
       call memevalf(n,nwinfeasx,ierr,memev)
       if ( ierr .ne. 0 ) return

       nwinfeasf = memev%f

       call updatesol(epsfeas,n,nwinfeasx,m,p,nwinfeasl,nwinfeasf,nwinfeascsupn,nwinfeasbdsvio,nwinfeasssupn, &
            nwinfeasnlpsupn,kktfound,kktf,kktx,kktlambda,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
            feasfound,feasf,feasx,feaslambda,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
            infeasfound,infeasf,infeasx,infeaslambda,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
            bestf,bestx,bestlambda,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn, &
            epscompl,epsopt,outiter,totiter,nwcalls,nwtotit)

       write(*,*) 'kkt:    ',kktfound,kktf,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn
       write(*,*) 'feas:   ',feasfound,feasf,feascsupn,feasbdsvio,feasssupn,feasnlpsupn
       write(*,*) 'infeas: ',infeasfound,infeasf,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn
       write(*,*) 'best:   ',bestf,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn
    end if
    
    call memevrestore(memev,memevb)

    write(*,*) 'Subproblem minimization.'
       
    outiter = outiter + 1

    if ( outiter .eq. 1 ) then
       epsopk = sqrt( epsopt )
       maxit = 10

       if ( rhoauto ) then
          call memevalf(n,x,ierr,memev)
          if ( ierr .ne. 0 ) return
          
          call memevalc(n,x,ierr,memev)
          if ( ierr .ne. 0 ) return
    
          sumc = 0.5d0 * ( sum( memev%scaledc(1:m) ** 2 ) + sum( max( 0.0d0, memev%scaledc(m+1:m+p) ) ** 2 ) )
          
          rho = 10.0d0 * max( 1.0d0, abs( memev%scaledf ) ) / max( 1.0d0, sumc )
       else
          rho = rhoini
       end if
       
       rho = max( rhoinimin, min( rho, rhoinimax ) )

    else
       epsopk = max( epsopt, 0.1d0 * epsopk )
       maxit = 5000

       if ( ( csupn .gt. epsfeas  .and. .not. csupn .le. rhofrac * csupnprev ) .or. &
            ( ssupn .gt. epscompl .and. .not. ssupn .le. rhofrac * ssupnprev ) ) then
          rho = rhomult * rho
       end if
    end if

    if ( rho .gt. 1.0d+20 ) then
       write(*,*) 'The penalty parameter is larger than 1e+20. Insistence is futile.'
       istop = 4
       return
    end if
    
    ftarget = - 1.0d+12
    hfixstr = .false.

    if ( memev%scale ) then
       lambda(1:m+p) = memev%sf * lambda(1:m+p) / memev%sc(1:m+p)
    end if

    if ( .not. all( lammin .le. lambda(1:m+p) .and. lambda(1:m+p) .le. lammax ) ) then
       lambda(1:m+p) = 0.0d0
    end if
    
    aldata%rho => rho
    aldata%lambda => lambda(1:m+p)
    aldata%workl => workl(1:m+p)
    aldata%workr => workr(1:m+p)
    aldata%memev => memev

    call gencanb(evalal,evalnal,evalhal,hfixstr,n,x,lind,lbnd,uind,ubnd, &
         al,nl,nlpsupn,lss,ftarget,epsopk,maxit,extallowed,inniter,ierr, &
         innistop,pdata=c_loc(aldata))
    
    if ( ierr .ne. 0 ) return
    
    totiter = totiter + inniter

    bdsvio = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    lambda(1:m+p) = lambda(1:m+p) + rho * memev%scaledc(1:m+p)

    lambda(m+1:m+p) = max( 0.0d0, lambda(m+1:m+p) )

    if ( memev%scale ) then
       lambda(1:m+p) = lambda(1:m+p) * memev%sc(1:m+p) / memev%sf
    end if

    csupnprev = csupn
    ssupnprev = ssupn
    
    csupn = max( 0.0d0, max( maxval( abs( memev%c(1:m) ) ), maxval( memev%c(m+1:m+p) ) ) )
    ssupn = max( 0.0d0, maxval( min( - memev%c(m+1:m+p), lambda(m+1:m+p) ) ) )

    write(*,*) 'outiter = ',outiter,' rho = ',rho,' csupn = ',csupn,' ssupn = ',ssupn,' nlpsupn = ',nlpsupn
    write(*,*) 'bdsvio = ',bdsvio,' (Should be zero here.)'
    
    if ( max( csupn, bdsvio ) .le. epsfeas ) then
       call memevalf(n,x,ierr,memev)
       if ( ierr .ne. 0 ) return

       f = memev%f

       call updatesol(epsfeas,n,x,m,p,lambda,f,csupn,bdsvio,ssupn,nlpsupn, &
            kktfound,kktf,kktx,kktlambda,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
            feasfound,feasf,feasx,feaslambda,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
            infeasfound,infeasf,infeasx,infeaslambda,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
            bestf,bestx,bestlambda,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn, &
            epscompl,epsopt,outiter,totiter,nwcalls,nwtotit)

       write(*,*) 'kkt:    ',kktfound,kktf,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn
       write(*,*) 'feas:   ',feasfound,feasf,feascsupn,feasbdsvio,feasssupn,feasnlpsupn
       write(*,*) 'infeas: ',infeasfound,infeasf,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn
       write(*,*) 'best:   ',bestf,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn
    
       if ( ssupn .le. epscompl .and. nlpsupn .le. epsopt ) then
          write(*,*) 'Success (after a subproblem minimization)!'
          istop = 0
          return
       end if
    end if
    
    if ( csupn .le. epsfeas .or. outiter .eq. 1 .or. csupn .lt. csupnprev ) then
       cnoprog = 0
    else
       nfeas(1:n) = 0.0d0

       workr(1:m) = memev%scaledc(1:m)
       workr(m+1:m+p) = max( 0.0d0, memev%scaledc(m+1:m+p) )
       
       workl(1:m+p) = workr(1:m+p) .ne. 0.0d0

       if ( any( workl(1:m+p) ) ) then
          call memevalj(n,x,m,p,workl,ierr,memev)
          if ( ierr .ne. 0 ) return
    
          do j = 1,m + p
             if ( workr(j) .ne. 0.0d0 ) then
                a = memev%jsta(j)
                b = a + memev%jlen(j) - 1
                nfeas(memev%jvar(a:b)) = nfeas(memev%jvar(a:b)) + workr(j) * memev%scaledjval(a:b)
             end if
          end do
       end if

       nfeas(1:n) = x(1:n) - nfeas(1:n)
       
       call project(n,lind,lbnd,uind,ubnd,nfeas)

       nfeas(1:n) = nfeas(1:n) - x(1:n)
       
       nfeaspsupn = maxval( abs( nfeas(1:n) ) )
       
       write(*,*) 'The required feasibility tolerance was not reached yet.'
       write(*,*) 'Moreover, in the current iteration there was not even simple deacrease in the feasibility.'
       write(*,*) 'This is the sup-norm of the continous projected gradient of 1/2 of the squared ', &
            '(scaled) infeasibility = ',nfeaspsupn

       if ( inniter .gt. 0 ) then
          cnoprog = cnoprog + 1
          if ( cnoprog .ge. 3 ) then
             write(*,*) 'Lack of progress in feasibility, zeroing the Langrange multipliers.'
             lambda(1:m+p) = 0.0d0
          end if
       end if
    end if
    
    if ( csupn .le. epsfeas ) then
       if ( innistop .eq. 0 ) then
          subfail = 0
       else
          subfail = subfail + 1
          if ( subfail .ge. 3 ) then
             write(*,*) 'Too many consecutive fails of the inner solver'
             istop = 2
             return
          end if
       end if
    end if

    go to 100

  end subroutine algencanb

  ! *****************************************************************
  ! *****************************************************************

  subroutine updatesol(epsfeas,n,x,m,p,lambda,f,csupn,bdsvio,ssupn,nlpsupn, &
       kktfound,kktf,kktx,kktlambda,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
       feasfound,feasf,feasx,feaslambda,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
       infeasfound,infeasf,infeasx,infeaslambda,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
       bestf,bestx,bestlambda,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn, &
       epscompl,epsopt,outiter,totiter,nwcalls,nwtotit)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(inout) :: kktfound,feasfound,infeasfound
    integer, intent(in) :: m,n,p,outiter,totiter,nwcalls,nwtotit
    real(kind=8), intent(in) :: epsfeas,f,csupn,bdsvio,ssupn,nlpsupn,epscompl,epsopt
    real(kind=8), intent(inout) :: kktf,kktcsupn,kktbdsvio,kktssupn,kktnlpsupn, &
         feasf,feascsupn,feasbdsvio,feasssupn,feasnlpsupn, &
         infeasf,infeascsupn,infeasbdsvio,infeasssupn,infeasnlpsupn, &
         bestf,bestcsupn,bestbdsvio,bestssupn,bestnlpsupn

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: lambda(m+p),x(n)
    real(kind=8), intent(inout) :: kktx(n),kktlambda(m+p),feasx(n),feaslambda(m+p),infeasx(n), &
         infeaslambda(m+p),bestx(n),bestlambda(m+p)

    ! PARAMETERS
    real(kind=8), parameter :: ftol = 1.0d-08

    if ( max( csupn, bdsvio ) .le. epsfeas ) then
       
       if ( ssupn .le. epscompl .and. nlpsupn .le. epsopt ) then

          if ( .not. kktfound .or. ( kktfound .and. .not. f .gt. kktf ) ) then
             kktfound = .true.
             
             kktf = f
                
             kktx(1:n) = x(1:n)
             kktlambda(1:m+p) = lambda(1:m+p)
                
             kktcsupn = csupn
             kktssupn = ssupn
             kktnlpsupn = nlpsupn
             kktbdsvio = bdsvio
          end if

       else

          if ( .not. feasfound .or. ( feasfound .and. .not. f .gt. feasf ) ) then
             feasfound = .true.
             
             feasf = f
                
             feasx(1:n) = x(1:n)
             feaslambda(1:m+p) = lambda(1:m+p)
                
             feascsupn = csupn
             feasssupn = ssupn
             feasnlpsupn = nlpsupn
             feasbdsvio = bdsvio
          end if

       end if

    else if ( .not. feasfound .and. ( .not. infeasfound .or. ( infeasfound .and. &
         max( csupn, bdsvio ) .le. max( infeascsupn, infeasbdsvio ) ) ) ) then
       infeasfound = .true.
             
       infeasf = f
                
       infeasx(1:n) = x(1:n)
       infeaslambda(1:m+p) = lambda(1:m+p)
                
       infeascsupn = csupn
       infeasssupn = ssupn
       infeasnlpsupn = nlpsupn
       infeasbdsvio = bdsvio
    end if
    
    if ( kktfound .and. feasfound ) then
       if ( isnan( feasf ) .or. ( .not. isnan( kktf ) .and. &
            .not. kktf .gt. feasf + ftol * max( 1.0d0, abs( feasf ) ) ) ) then
          bestf = kktf

          bestx(1:n) = kktx(1:n)
          bestlambda(1:m+p) = kktlambda(1:m+p)
                
          bestcsupn = kktcsupn
          bestssupn = kktssupn
          bestnlpsupn = kktnlpsupn
          bestbdsvio = kktbdsvio
       else
          bestf = feasf

          bestx(1:n) = feasx(1:n)
          bestlambda(1:m+p) = feaslambda(1:m+p)
          
          bestcsupn = feascsupn
          bestssupn = feasssupn
          bestnlpsupn = feasnlpsupn
          bestbdsvio = feasbdsvio
       end if
          
    else if ( kktfound ) then
       bestf = kktf

       bestx(1:n) = kktx(1:n)
       bestlambda(1:m+p) = kktlambda(1:m+p)
          
       bestcsupn = kktcsupn
       bestssupn = kktssupn
       bestnlpsupn = kktnlpsupn
       bestbdsvio = kktbdsvio

    else if ( feasfound ) then
       bestf = feasf
       
       bestx(1:n) = feasx(1:n)
       bestlambda(1:m+p) = feaslambda(1:m+p)
       
       bestcsupn = feascsupn
       bestssupn = feasssupn
       bestnlpsupn = feasnlpsupn
       bestbdsvio = feasbdsvio

    else
       bestf = infeasf
       
       bestx(1:n) = infeasx(1:n)
       bestlambda(1:m+p) = infeaslambda(1:m+p)
       
       bestcsupn = infeascsupn
       bestssupn = infeasssupn
       bestnlpsupn = infeasnlpsupn
       bestbdsvio = infeasbdsvio
    end if

    open(20,file='alsolver-interrupted-tabline.txt')
    write(20,9000) n,m+p,bestf,bestcsupn,bestssupn,bestnlpsupn,bestbdsvio,outiter,totiter,nwcalls,nwtotit
    close(20)

9000 format(2(1X,I6),1X,1P,D24.16,4(1X,1P,D7.1),4(1X,I8))

  end subroutine updatesol
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalal(n,x,al,ierr,aldataptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    real(kind=8), intent(out) :: al
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: aldataptr
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)

    ! LOCAL SCALARS
    integer :: m,p
    type(aldata_type), pointer :: aldata
    type(memev_type), pointer :: memev
    
    call c_f_pointer(aldataptr,aldata)
    memev => aldata%memev
    
    call memevalf(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    al = memev%scaledf

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    aldata%workr(1:m) = 0.5d0 * aldata%rho * ( memev%scaledc(1:m) + aldata%lambda(1:m) / aldata%rho ) ** 2
    aldata%workr(m+1:m+p) = 0.5d0 * aldata%rho * max( 0.0d0, memev%scaledc(m+1:m+p) + aldata%lambda(m+1:m+p) / aldata%rho ) ** 2

    al = al + sum( aldata%workr(1:m+p) )

  end subroutine evalal
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalnal(n,x,nal,ierr,aldataptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    integer, intent(inout) :: ierr
    type(c_ptr), optional, intent(in) :: aldataptr

    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    real(kind=8), intent(out) :: nal(n)

    ! LOCAL SCALARS
    integer :: a,b,j,m,p
    type(aldata_type), pointer :: aldata
    type(memev_type), pointer :: memev

    call c_f_pointer(aldataptr,aldata)
    memev => aldata%memev

    call memevalg(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    nal(1:n) = memev%scaledg(1:n)

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    aldata%workr(1:m+p) = aldata%lambda(1:m+p) + aldata%rho * memev%scaledc(1:m+p)
    aldata%workr(m+1:m+p) = max( 0.0d0, aldata%workr(m+1:m+p) )

    aldata%workl(1:m) = .true.
    aldata%workl(m+1:m+p) = aldata%workr(m+1:m+p) .gt. 0.0d0

    if ( any( aldata%workl(1:m+p) ) ) then
       call memevalj(n,x,m,p,aldata%workl,ierr,memev)
       if ( ierr .ne. 0 ) return
    
       do j = 1,m + p
          if ( aldata%workr(j) .ne. 0.0d0 ) then
             a = memev%jsta(j)
             b = a + memev%jlen(j) - 1
             nal(memev%jvar(a:b)) = nal(memev%jvar(a:b)) + aldata%workr(j) * memev%scaledjval(a:b)
          end if
       end do
    end if

  end subroutine evalnal
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine evalhal(n,x,lim,halnnz,halrow,halcol,halval,ierr,aldataptr)

    implicit none
    
    ! SCALAR ARGUMENTS
    integer, intent(in) :: lim,n
    integer, intent(inout) :: ierr
    integer, intent(out) :: halnnz
    type(c_ptr), optional, intent(in) :: aldataptr
    
    ! ARRAY ARGUMENTS
    real(kind=8), intent(in) :: x(n)
    integer, intent(out) :: halcol(lim),halrow(lim)
    real(kind=8), intent(out) :: halval(lim)

    ! LOCAL SCALARS
    integer :: a,b,dim,hlnnz,j,jdensit,jlenmax,k,l,m,nelem,p
    type(aldata_type), pointer :: aldata
    type(memev_type), pointer :: memev

    call c_f_pointer(aldataptr,aldata)
    memev => aldata%memev

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return
    
    m = memev%m
    p = memev%p

    aldata%workr(1:m+p) = aldata%lambda(1:m+p) + aldata%rho * memev%scaledc(1:m+p)
    aldata%workr(m+1:m+p) = max( 0.0d0, aldata%workr(m+1:m+p) )

    if ( memev%scale ) then
       aldata%workr(1:m+p) = memev%sc(1:m+p) * aldata%workr(1:m+p) / memev%sf
    end if
    
    if ( memev%pdatapresent ) then
       call memev%evalhl(n,x,m,p,aldata%workr,lim,.true.,hlnnz,halrow,halcol,halval,ierr,memev%pdata)
    else
       call memev%evalhl(n,x,m,p,aldata%workr,lim,.true.,hlnnz,halrow,halcol,halval,ierr)
    end if
          
    if ( memev%scale ) then
       halval(1:hlnnz) = memev%sf * halval(1:hlnnz)
    end if

    halnnz = hlnnz
    
    if ( ierr .ne. 0 ) return
          
    aldata%workl(1:m) = .true.
    aldata%workl(m+1:m+p) = aldata%workr(m+1:m+p) .gt. 0.0d0

    if ( any( aldata%workl(1:m+p) ) ) then
       call memevalj(n,x,m,p,aldata%workl,ierr,memev)
       if ( ierr .ne. 0 ) return

       jlenmax = maxval( memev%jlen(1:m+p), aldata%workl(1:m+p) )
       jdensit = sum( ( memev%jlen(1:m+p) / n ) * ( ( memev%jlen(1:m+p) + 1 ) / ( n + 1 ) ), aldata%workl(1:m+p) )
       
       if ( .not. aldata%avoidhal .and. ( n .le. aldata%nmin .or. &
            ( jlenmax .le. aldata%jdenmax * n .and. jdensit .le. aldata%hdenmax ) ) ) then

          halnnz = hlnnz
          
          outer: do j = 1,m + p
             if ( aldata%workl(j) ) then
                a = memev%jsta(j)
                b = a + memev%jlen(j) - 1
                do k = a,b
                   if ( memev%scaledjval(k) .ne. 0.0d0 ) then
                      if ( .not. memev%jsorted(j) ) then
                         do l = a,b
                            if ( memev%jvar(l) .ge. memev%jvar(k) .and. memev%scaledjval(l) .ne. 0.0d0 ) then
!!$                               if ( halnnz + 1 .gt. lim .or. .not. ( n .le. aldata%nmin .or. &
!!$                                    ( 2.0d0 * ( halnnz + 1 - hlnnz ) ) / n / ( n + 1 ) .le. aldata%hdenmax ) ) then
                               if ( halnnz + 1 .gt. lim ) then
                                  aldata%avoidhal = .true.
                                  exit outer
                               end if
                                  
                               halnnz = halnnz + 1
                               halval(halnnz) = aldata%rho * memev%scaledjval(k) * memev%scaledjval(l)
                               halcol(halnnz) = memev%jvar(l)
                               halrow(halnnz) = memev%jvar(k)
                            end if
                         end do
                      else
                         do l = k,b
                            if ( memev%scaledjval(l) .ne. 0.0d0 ) then
!!$                               if ( halnnz + 1 .gt. lim .or. .not. ( n .le. aldata%nmin .or. &
!!$                                    ( 2.0d0 * ( halnnz + 1 - hlnnz ) ) / n / ( n + 1 ) .le. aldata%hdenmax ) ) then
                               if ( halnnz + 1 .gt. lim ) then
                                  aldata%avoidhal = .true.
                                  exit outer
                               end if
                                  
                               halnnz = halnnz + 1
                               halval(halnnz) = aldata%rho * memev%scaledjval(k) * memev%scaledjval(l)
                               halcol(halnnz) = memev%jvar(l)
                               halrow(halnnz) = memev%jvar(k)
                            end if
                         end do
                      end if
                   end if
                end do
             end if
          end do outer

          if ( .not. aldata%avoidhal ) then
             call collapse(n,halnnz,halrow,halcol,halval)
             call removenull(halnnz,halrow,halcol,halval)
          end if
          
       end if

       if ( aldata%avoidhal .or. .not. ( n .le. aldata%nmin .or. &
            ( jlenmax .le. aldata%jdenmax * n .and. jdensit .le. aldata%hdenmax ) ) ) then

          dim = n
    
          halnnz = hlnnz
          
          do j = 1,m + p
             if ( aldata%workl(j) ) then
                nelem = memev%jlen(j)
             
                if ( halnnz + nelem + 1 .gt. lim ) then
                   write(*,*) 'In evalhal (algencan.f90): lack of memory!'
                   ierr = -98
                   return
                end if
             
                dim = dim + 1

                a = memev%jsta(j)
                b = a + nelem - 1

                halrow(halnnz+1:halnnz+nelem) = memev%jvar(a:b)
                halcol(halnnz+1:halnnz+nelem) = dim
                halval(halnnz+1:halnnz+nelem) = memev%scaledjval(a:b)

                halrow(halnnz+nelem+1) = dim
                halcol(halnnz+nelem+1) = dim
                halval(halnnz+nelem+1) = - 1.0d0 / aldata%rho

                halnnz = halnnz + nelem + 1
             end if
          end do
          
       end if
    end if

  end subroutine evalhal
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine collapse(n,halnnz,halrow,halcol,halval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: halnnz,n

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: halcol(halnnz),halrow(halnnz)
    real(kind=8), intent(inout) :: halval(halnnz)

    ! LOCAL SCALARS
    integer :: col,i,k,next,row
    
    ! LOCAL ARRAYS
    integer :: halnext(halnnz),halrowsta(n)
    real(kind=8) :: r(n)
    
    halrowsta(1:n) = 0
       
    do k = 1,halnnz
       row = halrow(k)
       halnext(k) = halrowsta(row)
       halrowsta(row) = k
    end do

    r(1:n) = 0.0d0
    
    do i = 1,n
       next = halrowsta(i)
       do while ( next .ne. 0 )
          col = halcol(next)
          r(col) = r(col) + halval(next)
          next = halnext(next)
       end do

       next = halrowsta(i)
       do while ( next .ne. 0 )
          col = halcol(next)
          halval(next) = r(col)
          r(col) = 0.0d0
          next = halnext(next)
       end do
    end do
       
  end subroutine collapse

  ! *****************************************************************
  ! *****************************************************************

  subroutine removenull(halnnz,halrow,halcol,halval)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(inout) :: halnnz

    ! ARRAY ARGUMENTS
    integer, intent(inout) :: halcol(halnnz),halrow(halnnz)
    real(kind=8), intent(inout) :: halval(halnnz)

    ! LOCAL SCALARS
    integer :: k,tmp
    
    ! LOCAL ARRAYS
    logical :: nonnull(halnnz)
    integer :: nnind(halnnz)

    nonnull(1:halnnz) = halval(1:halnnz) .ne. 0.0d0
    
    tmp = count( nonnull(1:halnnz) )
    nnind(1:tmp) = pack( (/ (k, k=1,halnnz) /), nonnull(1:halnnz) )
       
    halrow(1:tmp) = halrow(nnind(1:tmp))
    halcol(1:tmp) = halcol(nnind(1:tmp))
    halval(1:tmp) = halval(nnind(1:tmp))
    halnnz = tmp

  end subroutine removenull

end module bmalgencan
