module bmgencan

  use lss
  use iso_c_binding, only: c_ptr
  
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

     subroutine hsub(n,x,lim,hnnz,hrow,hcol,hval,ierr,pdata) bind(C)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n,lim
       real(kind=8), intent(in) :: x(n)
       integer, intent(out) :: hnnz
       integer, intent(out) :: hrow(lim),hcol(lim)
       real(kind=8), intent(out) :: hval(lim)
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine hsub

     subroutine stpsub(n,x,gsupn,inhdefstp,stp,ierr,pdata) bind(C)
       use iso_c_binding, only: c_ptr
       implicit none
       integer, intent(in) :: n
       real(kind=8), intent(in) :: gsupn
       real(kind=8), intent(in) :: x(n)
       logical, intent(out) :: inhdefstp,stp
       integer, intent(inout) :: ierr
       type(c_ptr), optional, intent(in) :: pdata
     end subroutine stpsub
  end interface

  private
  
  ! SCALAR PARAMETERS
  integer, parameter :: dimmax = 1000000, texmax = 20
  real(kind=8), parameter :: lrgstreal = huge( 1.0d0 ), macheps = epsilon( 1.0d0 ), &
       macheps12 = sqrt( macheps ), macheps14 = sqrt( macheps12 ), beta = 0.5d0, &
       gamma = 1.0d-04, eta = 1.0d+04, lspgmin = 1.0d-16, lspgmax = 1.0d+16, r = 0.1d0, &
       tau1 = 0.1d0, tau2 = 0.9d0, sgmsmall = 1.0d-08, sgmmin = 1.0d-08, sgmmax = 1.0d+16, &
       hmin = 1.0d-08, hmax = 1.0d+08

  public genunc, genuncb, gencan, gencanb, project

contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine genunc(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,f,g,gsupn, &
       ftarget,eps,maxit,extallowed,iter,ierr,istop,evalstp,pdata) bind(C, name="genunc")

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(hsub) :: evalh
    procedure(stpsub), optional :: evalstp

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,hfixstr
    integer, intent(in) :: hnnzmax,maxit,n
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: ftarget,eps
    real(kind=8), intent(out) :: f,gsupn
    type(c_ptr), optional, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    type(lss_type) :: lss

    call lssini(hnnzmax,lss)

    ierr = 0
    
    call genuncb(evalf,evalg,evalh,hfixstr,n,x,f,g,gsupn,lss,ftarget,eps, &
         maxit,extallowed,iter,ierr,istop,evalstp,pdata)

    call lssend(lss)
    
  end subroutine genunc

  ! *****************************************************************
  ! *****************************************************************

  subroutine genuncb(evalf,evalg,evalh,hfixstr,n,x,f,g,gsupn,lss,ftarget, &
       eps,maxit,extallowed,iter,ierr,istop,evalstp,pdata)

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(hsub) :: evalh
    procedure(stpsub), optional :: evalstp

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,hfixstr
    integer, intent(in) :: maxit,n
    integer, intent(inout) :: ierr
    integer, intent(out) :: istop,iter
    real(kind=8), intent(in) :: ftarget,eps
    real(kind=8), intent(out) :: f,gsupn
    type(c_ptr), optional, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: g(n)
    type(lss_type), intent(inout) :: lss

    ! istop =  0 means gsupn <= eps
    ! istop =  1 means gsupn <= eps^1/2 during 100 consecutive iterations
    ! istop =  2 means gsupn <= eps^1/4 during 5,000 consecutive iterations
    ! istop =  3 means gsupn <= eps^1/8 during 10,000 consecutive iterations
    ! istop =  6 means f <= ftarget
    ! istop =  8 means lack of progress (f did not improve in 3 consecuitove iterations)
    ! istop =  9 means iterates are diverging ( || xk ||_2 >= huge( 1.0d0 ) )
    ! istop = 10 means iter >= maxit
    ! istop = 20 means user-provided stopping criterion satisfied
    
    ! LOCAL SCALARS
    logical :: extinhibited,extrecom,inhdefstp,lssana,rankdef,sgminiundef,ustp
    integer :: dim,extfails,i,k,navoidcy,nfnoprogr,nneigv,tex
    real(kind=8) :: eps12,eps14,eps18,fbest,fbtrial,fext,fref,ftrial, &
         gtd,gsupnb,sgmini,sgmtrial,t,ttmp
  
    ! LOCAL ARRAYS
    integer :: countg(4),idiag(n)
    real(kind=8) :: d(n),diagv(n),gtrial(n),origv(n),rhs(dimmax),xbtrial(n), &
         xext(n),xref(n),xtrial(n)
  
    ! ==========================================================================
    ! Initialization
    ! ==========================================================================

    iter = 0

    extfails = 0
    extinhibited = .false.
    
    navoidcy = 0
    nfnoprogr = 0
    fbest = lrgstreal
    gsupnb = lrgstreal

    countg(1:4) = 0
    
    eps12 = sqrt( eps )
    eps14 = sqrt( eps12 )
    eps18 = sqrt( eps14 )

    sgminiundef = .true.
    
    call evalf(n,x,f,ierr,pdata)
    if ( ierr .ne. 0 ) return

    call evalg(n,x,g,ierr,pdata)
    if ( ierr .ne. 0 ) return

    ! ==========================================================================
    ! Main loop
    ! ==========================================================================
  
100 continue
  
    ! Check stopping criterion

    gsupn = maxval( abs( g(1:n) ) )

    write(*,*) 'iter = ',iter,' f = ',f,' gsupn = ',gsupn
  
    open(20,file='solver-interrupted-tabline.txt')
    write(20,9000) n,f,gsupn,iter
    close(20)

    inhdefstp = .false.
    
    if ( present( evalstp ) ) then
       call evalstp(n,x,gsupn,inhdefstp,ustp,ierr,pdata)
       if ( ierr .ne. 0 ) return

       if ( ustp ) then
          write(*,*) 'User-defined stopping criterion satisfied'
          istop = 20
          return
       end if
    end if
    
    if ( gsupn .le. eps ) then
       countg(1:4) = countg(1:4) + 1
    else
       countg(1:1) = 0
       if ( gsupn .le. eps12 ) then
          countg(2:4) = countg(2:4) + 1
       else
          countg(2:2) = 0
          if ( gsupn .le. eps14 ) then
             countg(3:4) = countg(3:4) + 1
          else
             countg(3:3) = 0
             if ( gsupn .le. eps18 ) then
                countg(4:4) = countg(4:4) + 1
             else
                countg(4:4) = 0
             end if
          end if
       end if
    end if

    if( .not. inhdefstp )  then
       if ( countg(1) .ge. 1 ) then
          write(*,*) 'Stopping criterion gsupn <= eps satisfied'
          istop = 0
          return
       else if ( countg(2) .ge. 100 ) then
          write(*,*) 'Stopping because gsupn <= eps^1/2 during 100 consecutive iterations'
          istop = 1
          return
       else if ( countg(3) .ge. 5000 ) then
          write(*,*) 'Stopping because gsupn <= eps^1/4 during 5,000 consecutive iterations'
          istop = 2
          return
       else if ( countg(4) .ge. 10000 ) then
          write(*,*) 'Stopping because gsupn <= eps^1/8 during 10,000 consecutive iterations'
          istop = 3
          return
       end if
  
       if ( f .le. ftarget ) then
          write(*,*) 'Stopping criterion f <= ftarget satisfied'
          istop = 6
          return
       end if

       if ( f .lt. fbest ) then
          nfnoprogr = 0
       
       else
          nfnoprogr = nfnoprogr + 1
          if ( nfnoprogr .gt. 3 ) then
             write(*,*) 'Lack of progress.'
             istop = 8
             return
          end if
       end if
       
       if ( f .lt. fbest .or. ( f .eq. fbest .and. gsupn .lt. gsupnb ) ) then
          fbest = f
          gsupnb = gsupn
       end if

       if ( norm2( x(1:n) ) .ge. lrgstreal ) then
          istop = 9
          return
       end if
  
       if ( iter .ge. maxit ) then
          istop = 10
          return
       end if
    end if

    ! iterate

    iter = iter + 1

    extrecom = .false.
    
    fbtrial = lrgstreal

    ! Compute Hessian

    call evalh(n,x,lss%nnzmax,lss%matrix%ne,lss%matrix%row,lss%matrix%col,lss%matrix%val,ierr,pdata)
    if ( ierr .ne. 0 ) return

    dim = max( n, maxval(lss%matrix%row(1:lss%matrix%ne)), maxval(lss%matrix%col(1:lss%matrix%ne)) )
      
    if ( dim .gt. dimmax )  then
       write(*,*) 'Increase dimmax in bmgencan to at least ',dim,' and re-run.'
       stop
    end if

!!$    lss%matrix%n = n

    lss%matrix%n = dim

    ! Set pointers to the diagonal elements

    idiag(1:n) = 0
      
    do k = 1,lss%matrix%ne
       if ( lss%matrix%row(k) .eq. lss%matrix%col(k) .and. lss%matrix%row(k) .le. n ) then
          if ( idiag(lss%matrix%row(k)) .eq. 0 ) then
             idiag(lss%matrix%row(k)) = k
             origv(lss%matrix%row(k)) = lss%matrix%val(k)
             diagv(lss%matrix%row(k)) = lss%matrix%val(k)
          else
             diagv(lss%matrix%row(k)) = diagv(lss%matrix%row(k)) + lss%matrix%val(k)
          end if
       end if
    end do
       
    do i = 1,n
       if ( idiag(i) .eq. 0 ) then
          if ( lss%matrix%ne + 1 .gt. lss%nnzmax ) then
             write(*,*) 'In genuncb (gencan.f90): lack of memory!'
             ierr = -98
             return
          end if
             
          lss%matrix%ne = lss%matrix%ne + 1
             
          lss%matrix%row(lss%matrix%ne) = i
          lss%matrix%col(lss%matrix%ne) = i
          lss%matrix%val(lss%matrix%ne) = 0.0d0
             
          idiag(i) = lss%matrix%ne
          origv(i) = 0.0d0
          diagv(i) = 0.0d0
       end if
    end do

    ! Compute search direction
    
    ! Try Newton direction first (i.e. unmodified Hessian)

    sgmtrial = 0.0d0

    lssana = iter .eq. 1 .or. .not. hfixstr
    call lssfac(lss,lssana,nneigv,rankdef)
    
    if ( nneigv .ne. dim - n .or. rankdef ) then
       write(*,*) 'The Hessian has negative eigenvalues or it is rank defficient.'
       go to 200
    end if
    
    rhs(1:n) = - g(1:n)
    rhs(n+1:dim) = 0.0d0

    call lsssol(lss,dim,rhs)

    d(1:n) = rhs(1:n)
    
    write(*,*) 'Newton direction dsupn = ',maxval( abs( d(1:n) ) ),' deucn = ',norm2( d(1:n) )

    go to 300

    ! Try with H + sigma Id until obtaining desired inertia
    
200 continue
             
    if ( sgminiundef ) then
       sgminiundef = .false.
       
       sgmini = sgmsmall * max( hmin, min( maxval( abs( diagv(1:n) ) ), hmax ) )
       sgmini = max( sgmmin, min( sgmini, sgmmax ) )
    end if
    
    sgmtrial = sgmini
       
210 continue

    if ( .not. sgmtrial .le. lrgstreal ) then
       istop = 11
       return
    end if
    
    write(*,*) 'We try now with Hessian + sgmtrial Id, sgmtrial = ',sgmtrial
    
    lss%matrix%val(idiag(1:n)) = origv(1:n) + sgmtrial
       
    lssana = .false.
    call lssfac(lss,lssana,nneigv,rankdef)

    if ( nneigv .ne. dim - n .or. rankdef ) then
       write(*,*) 'The Hessian + sigma Id has negative eigenvalues or it is rank defficient; so we increase sigma.'
       if ( 10.0d0 * sgmtrial .le. lrgstreal ) then
          sgmtrial = 10.0d0 * sgmtrial
          go to 210
       else
          write(*,*) 'sigma is already too big.'
       end if
    end if
    
    rhs(1:n) = - g(1:n)
    rhs(n+1:dim) = 0.0d0

    call lsssol(lss,dim,rhs)

    d(1:n) = rhs(1:n)

    write(*,*) 'dsupn = ',maxval( abs( d(1:n) ) ),' deucn = ',norm2( d(1:n) )
      
    if ( norm2( d(1:n) ) .gt. eta * max( 1.0d0, norm2( x(1:n) ) ) ) then
       write(*,*) 'The norm of the direction is too big; so we skip it and increase sigma.'
       if ( 10.0d0 * sgmtrial .le. lrgstreal ) then
          sgmtrial = 10.0d0 * sgmtrial
          go to 210
       else
          write(*,*) 'sigma is already too big.'
       end if
    end if
          
    ! Backtracking

300 continue
    
    gtd = dot_product( g(1:n), d(1:n) )

    write(*,*) 'gtd = ',gtd
      
    t = 1.0d0

310 continue
      
    xtrial(1:n) = x(1:n) + t * d(1:n)

    call evalf(n,xtrial,ftrial,ierr,pdata)
    if ( ierr .ne. 0 ) return
    
    write(*,*) 'in backtracking, t = ',t,' ftrial = ',ftrial

    if ( ftrial .lt. fbtrial ) then
       fbtrial = ftrial
       xbtrial(1:n) = xtrial(1:n)
    end if
      
    if ( t .eq. 0.0d0 ) then
       write(*,*) 'null step in backtracking'
       xtrial(1:n) = x(1:n)
       ftrial = f
       go to 400
       
    else if ( ftrial .le. ftarget ) then
       write(*,*) 'ftrial <= ftarget'
       go to 400
       
    else if ( ftrial .le. f + t * gamma * gtd ) then
       write(*,*) 'Armijo holds'
       if ( t .eq. 1.0d0 ) then
          extrecom = .true.
       end if
       go to 400

    else if ( t .eq. 1.0d0 .and. sgmtrial .eq. 0.0d0 .and. &
         ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. &
         norm2( d(1:n) ) .le. macheps14 * norm2( x(1:n) ) ) then
       if ( navoidcy .lt. 2 ) then
          navoidcy = navoidcy + 1
          write(*,*) 'The Newton step did NOT give sufficient decrease.'
          write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
          write(*,*) 'Therefore, we accept it anyway.'
          go to 400
       end if
    end if

    ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
      
    if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
       t = 0.5d0 * t
    else
       t = ttmp
    end if
    
    go to 310

400 continue
    
    if ( fbtrial .lt. ftrial ) then
       ftrial = fbtrial
       xtrial(1:n) = xbtrial(1:n)
    end if
       
    if ( sgmtrial .eq. 0.0d0 ) then
       sgmini = min( max( sgmmin, 0.5d0 * sgmini   ), sgmmax )
    else
       sgmini = min( max( sgmmin, 0.5d0 * sgmtrial ), sgmmax )
    end if
    
    call evalg(n,xtrial,gtrial,ierr,pdata)
    if ( ierr .ne. 0 ) return

    ! Extrapolate if recommended
    
    if ( extallowed .and. .not. extinhibited .and. extrecom ) then
       if ( .not. maxval( abs( gtrial(1:n) ) ) .le. eps ) then
          if ( dot_product( gtrial(1:n), d(1:n) ) .le. beta * gtd ) then
             fref = ftrial
             xref(1:n) = xtrial(1:n)
    
             tex = 1
410          if ( tex .le. texmax ) then
                if ( .not. fref .le. ftarget ) then
                   xext(1:n) = x(1:n) + 2.0d0 ** tex * ( xtrial(1:n) - x(1:n) )
                   
                   call evalf(n,xext,fext,ierr,pdata)
                   if ( ierr .ne. 0 ) return
                   
                   if ( fext .lt. fref ) then
                      fref = fext
                      xref(1:n) = xext(1:n)
                      
                      write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' gave descent fext = ',fext
                      
                      tex = tex + 1
                      go to 410
                   else
                      write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' failed.'
                      if ( tex .gt. 1 )  then
                         extfails = 0
                      else
                         extfails = extfails + 1
                         if ( extfails .eq. 3 ) then
                            extinhibited = .true.
                         end if
                      end if
                   end if
                else
                   write(*,*) 'fref <= ftarget'
                end if
             else
                write(*,*) 'tex > texmax (tex =',tex,', texmax = ',texmax,')'
             end if

             if ( fref .lt. ftrial ) then
                ftrial = fref
                xtrial(1:n) = xref(1:n)
                
                call evalg(n,xtrial,gtrial,ierr,pdata)
                if ( ierr .ne. 0 ) return
             end if
          end if
       end if
    end if

    f = ftrial
    x(1:n) = xtrial(1:n)
    g(1:n) = gtrial(1:n)

    go to 100
  
9000 format(1X,I6,1X,1P,D24.16,1X,1P,D7.1,1X,I8)

  end subroutine genuncb
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine gencan(evalf,evalg,evalh,hnnzmax,hfixstr,n,x,lind, &
       lbnd,uind,ubnd,f,g,gpsupn,ftarget,eps,maxit,extallowed,iter, &
       ierr,istop,evalstp,pdata) bind(C, name="gencan")

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(hsub) :: evalh
    procedure(stpsub), optional :: evalstp

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,hfixstr
    integer, intent(in) :: hnnzmax,maxit,n
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: ftarget,eps
    real(kind=8), intent(out) :: f,gpsupn
    type(c_ptr), optional, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: g(n)

    ! LOCAL SCALARS
    type(lss_type) :: lss

    call lssini(hnnzmax,lss)

    ierr = 0
    
    call gencanb(evalf,evalg,evalh,hfixstr,n,x,lind,lbnd,uind,ubnd, &
         f,g,gpsupn,lss,ftarget,eps,maxit,extallowed,iter,ierr,istop, &
         evalstp,pdata)

    call lssend(lss)
    
  end subroutine gencan

  ! *****************************************************************
  ! *****************************************************************

  subroutine gencanb(evalf,evalg,evalh,hfixstr,n,x,lind,lbnd,uind,ubnd, &
       f,g,gpsupn,lss,ftarget,eps,maxit,extallowed,iter,ierr,istop,evalstp,pdata)

    implicit none

    ! PROCEDURE ARGUMENTS
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(hsub) :: evalh
    procedure(stpsub), optional :: evalstp

    ! SCALAR ARGUMENTS
    logical, intent(in) :: extallowed,hfixstr
    integer, intent(in) :: maxit,n
    integer, intent(inout) :: ierr
    integer, intent(out) :: istop,iter
    real(kind=8), intent(in) :: ftarget,eps
    real(kind=8), intent(out) :: f,gpsupn
    type(c_ptr), optional, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)
    real(kind=8), intent(out) :: g(n)
    type(lss_type), intent(inout) :: lss
  
    ! Solves the box-constrained minimization problem
    !
    !        Minimize f(x) subject to l <= x <= u
    !
    ! using the method described in
    !
    ! E. G. Birgin and J. M. Martinez, ''Large-scale active-set box-
    ! constrained optimization method with spectral projected
    ! gradients'', Computational Optimization and Applications 23, pp.
    ! 101-125, 2002.

    ! istop =  0 means gsupn <= eps
    ! istop =  1 means gsupn <= eps^1/2 during 100 consecutive iterations
    ! istop =  2 means gsupn <= eps^1/4 during 5,000 consecutive iterations
    ! istop =  3 means gsupn <= eps^1/8 during 10,000 consecutive iterations
    ! istop =  6 means f <= ftarget
    ! istop =  8 means lack of progress (f did not improve in 3 consecuitove iterations)
    ! istop =  9 means iterates are diverging ( || xk ||_2 >= huge( 1.0d0 ) )
    ! istop = 10 means iter >= maxit
    ! istop = 11 means sigma trial is too large (this should never happen)
    ! istop = 20 means user-provided stopping criterion satisfied
    
    ! LOCAL SCALARS
    logical :: extinhibited,extrecom,infrontier,inhdefstp,newface,sgminiundef,rankdef,ustp
    integer :: extfails,i,navoidcy,nfnoprogr,nfr,nneigv,tex
    real(kind=8) :: eps12,eps14,eps18,fbest,fbtrial,fext,fref,ftrial,gisupn,gtd, &
         gpsupnb,lspg,sts,sty,sgmini,sgmtrial,t,ttmp
  
    ! LOCAL ARRAYS
    logical :: frvar(n),fvpre(n)
    integer :: countg(4),ifr(n)
    real(kind=8) :: d(n),gp(n),gprev(n),gptrial(n),gtrial(n),s(n),xbtrial(n), &
         xext(n),xprev(n),xref(n),xtrial(n)
  
    ! ==========================================================================
    ! Initialization
    ! ==========================================================================

    iter = 0

    extfails = 0
    extinhibited = .false.
    
    navoidcy = 0
    nfnoprogr = 0
    fbest = lrgstreal
    gpsupnb = lrgstreal

    countg(1:4) = 0
    
    eps12 = sqrt( eps )
    eps14 = sqrt( eps12 )
    eps18 = sqrt( eps14 )

    sgminiundef = .true.

    call project(n,lind,lbnd,uind,ubnd,x)

    call evalf(n,x,f,ierr,pdata)
    if ( ierr .ne. 0 ) return

    call evalg(n,x,g,ierr,pdata)
    if ( ierr .ne. 0 ) return

    ! ==========================================================================
    ! Main loop
    ! ==========================================================================
  
100 continue
  
    ! Check stopping criterion

    gp(1:n) = x(1:n) - g(1:n)

    call project(n,lind,lbnd,uind,ubnd,gp)

    gp(1:n) = gp(1:n) - x(1:n)
  
    gpsupn = maxval( abs( gp(1:n) ) )

    write(*,*) 'iter = ',iter,' f = ',f,' gpsupn = ',gpsupn
  
    open(20,file='solver-interrupted-tabline.txt')
    write(20,9000) n,f,gpsupn,iter
    close(20)

    inhdefstp = .false.
    
    if ( present( evalstp ) ) then
       call evalstp(n,x,gpsupn,inhdefstp,ustp,ierr,pdata)
       if ( ierr .ne. 0 ) return

       if ( ustp ) then
          write(*,*) 'User-defined stopping criterion satisfied'
          istop = 20
          return
       end if
    end if
    
    if ( gpsupn .le. eps ) then
       countg(1:4) = countg(1:4) + 1
    else
       countg(1:1) = 0
       if ( gpsupn .le. eps12 ) then
          countg(2:4) = countg(2:4) + 1
       else
          countg(2:2) = 0
          if ( gpsupn .le. eps14 ) then
             countg(3:4) = countg(3:4) + 1
          else
             countg(3:3) = 0
             if ( gpsupn .le. eps18 ) then
                countg(4:4) = countg(4:4) + 1
             else
                countg(4:4) = 0
             end if
          end if
       end if
    end if

    if ( .not. inhdefstp ) then
       if ( countg(1) .ge. 1 ) then
          write(*,*) 'Stopping criterion gpsupn <= eps satisfied'
          istop = 0
          return
       else if ( countg(2) .ge. 100 ) then
          write(*,*) 'Stopping because gpsupn <= eps^1/2 during 100 consecutive iterations'
          istop = 1
          return
       else if ( countg(3) .ge. 5000 ) then
          write(*,*) 'Stopping because gpsupn <= eps^1/4 during 5,000 consecutive iterations'
          istop = 2
          return
       else if ( countg(4) .ge. 10000 ) then
          write(*,*) 'Stopping because gpsupn <= eps^1/8 during 10,000 consecutive iterations'
          istop = 3
          return
       end if
       
       if ( f .le. ftarget ) then
          write(*,*) 'Stopping criterion f <= ftarget satisfied'
          istop = 6
          return
       end if

       if ( f .lt. fbest ) then
          nfnoprogr = 0
       
       else
          nfnoprogr = nfnoprogr + 1
          if ( nfnoprogr .gt. 3 ) then
             write(*,*) 'Lack of progress.'
             istop = 8
             return
          end if
       end if
       
       if ( f .lt. fbest .or. ( f .eq. fbest .and. gpsupn .lt. gpsupnb ) ) then
          fbest = f
          gpsupnb = gpsupn
       end if

       if ( norm2( x(1:n) ) .ge. lrgstreal ) then
          write(*,*) 'x components are too large'
          istop = 9
          return
       end if
  
       if ( iter .ge. maxit ) then
          write(*,*) 'Maximum of iterations reached'
          istop = 10
          return
       end if
    end if
    
    ! opt by an inner-to-the-face or a leaving-face iteration

    frvar(1:n) = ( .not. lind(1:n) .or. ( x(1:n) - lbnd(1:n) ) / max( 1.0d0, lbnd(1:n) ) .gt. 1.0d-12 ) .and. &
                 ( .not. uind(1:n) .or. ( ubnd(1:n) - x(1:n) ) / max( 1.0d0, ubnd(1:n) ) .gt. 1.0d-12 )

    newface = .false.
    if ( iter .eq. 0 .or. any( frvar(1:n) .neqv. fvpre(1:n) .or. &
         ( .not. frvar(1:n) .and. x(1:n) .ne. xprev(1:n) ) ) ) then
       newface = .true.
    end if

    write(*,*) 'Value of newface = ',newface
  
    ! iterate

    iter = iter + 1

    extrecom = .false.
    
    xprev(1:n) = x(1:n)
    gprev(1:n) = g(1:n)
    fvpre(1:n) = frvar(1:n)
  
    gisupn = max( 0.0d0, maxval( abs( gp(1:n) ), frvar(1:n) ) )

    write(*,*) 'gisupn = ',gisupn

    nfr = count( frvar(1:n) )

    write(*,*) 'Number of free variables: ',nfr
  
    if ( gisupn .ge. r * gpsupn ) then

       ifr(1:nfr) = pack( (/ (i, i=1,n) /), frvar(1:n) )

       write(*,*) 'An inner-to-the-face iteration will be done'

       call newtonls()
       
       if ( ierr .ne. 0 ) return

       if ( .not. sgmtrial .le. lrgstreal ) then
          write(*,*) 'sgmtrial is too large'
          istop = 11
          return
       end if
    
       if ( sgmtrial .eq. 0.0d0 ) then
          sgmini = min( max( sgmmin, 0.5d0 * sgmini   ), sgmmax )
       else
          sgmini = min( max( sgmmin, 0.5d0 * sgmtrial ), sgmmax )
       end if
    
       if ( fbtrial .lt. ftrial ) then
          ftrial = fbtrial
          xtrial(1:n) = xbtrial(1:n)
       end if
       
    else

       write(*,*) 'A leaving-face iteration will be done'

       if ( iter .eq. 1 .or. sty .le. 0.0d0 ) then
          lspg = max( 1.0d0, norm2( x(1:n) ) ) / norm2( gp(1:n) )
       else
          lspg = sts / sty
       end if

       lspg = max( lspgmin, min( lspg, lspgmax ) )

       t = 1.0d0

       xtrial(1:n) = x(1:n) - lspg * g(1:n)
       call project(n,lind,lbnd,uind,ubnd,xtrial)

       d(1:n) = xtrial(1:n) - x(1:n)

       gtd = dot_product( g(1:n), d(1:n) )

200    continue
     
       call evalf(n,xtrial,ftrial,ierr,pdata)
       if ( ierr .ne. 0 ) return

       if ( t .eq. 0.0d0 ) then
          write(*,*) 'null step in backtracking'
          xtrial(1:n) = x(1:n)
          ftrial = f
          
       else if ( ftrial .le. ftarget ) then
          write(*,*) 'ftrial <= ftarget'

       else if ( .not. ( ftrial .le. f + t * gamma * gtd ) ) then
          ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
          
          if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
             t = 0.5d0 * t
          else
             t = ttmp
          end if

          xtrial(1:n) = x(1:n) + t * d(1:n)

          go to 200
       end if

    end if

    call evalg(n,xtrial,gtrial,ierr,pdata)
    if ( ierr .ne. 0 ) return

    if ( extallowed .and. .not. extinhibited .and. extrecom ) then
       gptrial(1:n) = xtrial(1:n) - gtrial(1:n)
       call project(n,lind,lbnd,uind,ubnd,gptrial)
       gptrial(1:n) = gptrial(1:n) - xtrial(1:n)
  
       if ( .not. maxval( abs( gptrial(1:n) ) ) .le. eps ) then
          if ( infrontier .or. dot_product( gtrial(ifr(1:nfr)), d(1:nfr) ) .le. beta * gtd ) then
             fref = ftrial
             xref(1:n) = xtrial(1:n)
    
             tex = 1
310          if ( tex .le. texmax ) then
                if ( .not. fref .le. ftarget ) then
                   xext(1:n) = x(1:n) + 2.0d0 ** tex * ( xtrial(1:n) - x(1:n) )

                   call project(n,lind,lbnd,uind,ubnd,xext)

                   call evalf(n,xext,fext,ierr,pdata)
                   if ( ierr .ne. 0 ) return

                   if ( fext .lt. fref ) then
                      fref = fext
                      xref(1:n) = xext(1:n)
                      
                      write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' gave descent fext = ',fext
                      
                      tex = tex + 1
                      go to 310
                   else
                      write(*,*) 'Extrapolation with 2^tex = ',2.0d0 ** tex,' failed.'
                      if ( tex .gt. 1 )  then
                         extfails = 0
                      else
                         extfails = extfails + 1
                         if ( extfails .eq. 3 ) then
                            extinhibited = .true.
                         end if
                      end if
                   end if
                else
                   write(*,*) 'fref <= ftarget'
                end if
             else
                write(*,*) 'tex > texmax (tex =',tex,', texmax = ',texmax,')'
             end if

             if ( fref .lt. ftrial ) then
                ftrial = fref
                xtrial(1:n) = xref(1:n)
                
                call evalg(n,xtrial,gtrial,ierr,pdata)
                if ( ierr .ne. 0 ) return
             end if
          end if
       end if
    end if

    f = ftrial
    x(1:n) = xtrial(1:n)
    g(1:n) = gtrial(1:n)
    
    s(1:n) = x(1:n) - xprev(1:n)
    sty = dot_product( s(1:n), g(1:n) - gprev(1:n) )
    sts = sum( s(1:n) ** 2 )
  
    go to 100
  
9000 format(1X,I6,1X,1P,D24.16,1X,1P,D7.1,1X,I8)

  contains
    
    ! *****************************************************************
    ! *****************************************************************
    
    subroutine newtonls()

      implicit none
      
      ! LOCAL SCALARS
      logical :: lssana
      integer :: dim,i,k,nkeepel,rbdnnz
      real(kind=8) :: amax
      
      ! LOCAL ARRAYS
      integer :: idiag(nfr),rbdind(nfr),renumbered(dimmax)
      real(kind=8) :: amaxv(nfr),origv(nfr),rhs(dimmax),diagv(nfr)
      character :: amaxt(nfr),rbdtype(nfr)

      fbtrial = lrgstreal
      
      ! Compute Hessian (in the Martínez-Lucio format if the problem
      ! corresponds to the subproblem of an Augmented Lagrangian
      ! subproblem)
      
      call evalh(n,x,lss%nnzmax,lss%matrix%ne,lss%matrix%row,lss%matrix%col,lss%matrix%val,ierr,pdata)
      if ( ierr .ne. 0 ) return

      ! Shrink Hessian eliminating rows and columns related to fixed
      ! variables (dimension must be redefined and remaining elements
      ! renumbered)

      dim = max( n, maxval(lss%matrix%row(1:lss%matrix%ne)), maxval(lss%matrix%col(1:lss%matrix%ne)) )
      
      if ( dim .gt. dimmax )  then
         write(*,*) 'Increase dimmax in bmgencan to at least ',dim,' and re-run.'
         stop
      end if

      lss%matrix%n = dim - n + nfr

      renumbered(1:n) = unpack( (/ (i,i=1,nfr) /), frvar(1:n), 0 )
      renumbered(n+1:dim) = (/ (nfr+i,i=1,dim-n) /)

      where ( renumbered(lss%matrix%row(1:lss%matrix%ne)) .ne. 0 .and. &
              renumbered(lss%matrix%col(1:lss%matrix%ne)) .ne. 0 )
         lss%matrix%row(1:lss%matrix%ne) = renumbered(lss%matrix%row(1:lss%matrix%ne))
         lss%matrix%col(1:lss%matrix%ne) = renumbered(lss%matrix%col(1:lss%matrix%ne))
      elsewhere
         lss%matrix%row(1:lss%matrix%ne) = 0
         lss%matrix%col(1:lss%matrix%ne) = 0
      end where

      nkeepel = count( lss%matrix%row(1:lss%matrix%ne) .ne. 0 )

      lss%matrix%val(1:nkeepel) = pack( lss%matrix%val(1:lss%matrix%ne), lss%matrix%row(1:lss%matrix%ne) .ne. 0 )
      
      lss%matrix%row(1:nkeepel) = pack( lss%matrix%row(1:lss%matrix%ne), lss%matrix%row(1:lss%matrix%ne) .ne. 0 )
      lss%matrix%col(1:nkeepel) = pack( lss%matrix%col(1:lss%matrix%ne), lss%matrix%col(1:lss%matrix%ne) .ne. 0 )

      lss%matrix%ne = nkeepel
      
      ! Set pointers to the diagonal elements

      idiag(1:nfr) = 0
      
      do k = 1,lss%matrix%ne
         if ( lss%matrix%row(k) .eq. lss%matrix%col(k) .and. lss%matrix%row(k) .le. nfr ) then
            if ( idiag(lss%matrix%row(k)) .eq. 0 ) then
               idiag(lss%matrix%row(k)) = k
               origv(lss%matrix%row(k)) = lss%matrix%val(k)
               diagv(lss%matrix%row(k)) = lss%matrix%val(k)
            else
               diagv(lss%matrix%row(k)) = diagv(lss%matrix%row(k)) + lss%matrix%val(k)
            end if
         end if
      end do

      do i = 1,nfr
         if ( idiag(i) .eq. 0 ) then
            if ( lss%matrix%ne + 1 .gt. lss%nnzmax ) then
               write(*,*) 'In newtonls (gencan.f90): lack of memory!'
               ierr = -98
               return
            end if
          
            lss%matrix%ne = lss%matrix%ne + 1
          
            lss%matrix%row(lss%matrix%ne) = i
            lss%matrix%col(lss%matrix%ne) = i
            lss%matrix%val(lss%matrix%ne) = 0.0d0
            
            idiag(i) = lss%matrix%ne
            origv(i) = 0.0d0
            diagv(i) = 0.0d0
         end if
      end do

      ! Compute search direction
    
      ! Try Newton direction first (i.e. unmodified Hessian)

      sgmtrial = 0.0d0

      lssana = newface .or. .not. hfixstr

      call lssfac(lss,lssana,nneigv,rankdef)

      write(*,*) 'We are expecting nfr = ',nfr,' positive eigenvalues and lss%matrix%n - nfr = ', &
           lss%matrix%n - nfr,' negative eigenvalues.'
      
      write(*,*) 'In the unmodified Newtonian matrix, we got nneigv = ',nneigv,' rankdef = ',rankdef
      
      if ( nneigv .ne. lss%matrix%n - nfr .or. rankdef ) then
         go to 200
      end if

      write(*,*) 'We are done with the matrix inertia.'
      
      rhs(1:nfr) = - g(ifr(1:nfr))
      rhs(nfr+1:lss%matrix%n) = 0.0d0
      
      call lsssol(lss,lss%matrix%n,rhs)

      d(1:nfr) = rhs(1:nfr)

      write(*,*) 'Newton direction dsupn = ',maxval( abs( d(1:nfr) ) ),' deucn = ',norm2( d(1:nfr) )

      go to 300
      
      ! Try with H + sigma Id until obtaining desired inertia
    
200   continue
      
      if ( sgminiundef ) then
         sgminiundef = .false.
         
         sgmini = sgmsmall * max( hmin, min( maxval( abs( diagv(1:nfr) ) ), hmax ) )
         sgmini = max( sgmmin, min( sgmini, sgmmax ) )
      end if
      
      sgmtrial = sgmini

210   continue

      if ( .not. sgmtrial .le. lrgstreal ) then
         return
      end if

      write(*,*) 'We will try adding sgmtrial = ',sgmtrial,' to the NW.'
      
      lss%matrix%val(idiag(1:nfr)) = origv(1:nfr) + sgmtrial

      lssana = .false.
      
      call lssfac(lss,lssana,nneigv,rankdef)

      write(*,*) 'We got, nneigv = ',nneigv,' rankdef = ',rankdef
      
      if ( nneigv .ne. lss%matrix%n - nfr .or. rankdef ) then
         if ( 10.0d0 * sgmtrial .le. lrgstreal ) then
            sgmtrial = 10.0d0 * sgmtrial
            go to 210
         else
            write(*,*) 'sigma is already too big.'
         end if
      end if

      write(*,*) 'We are done with the matrix inertia.'
      
      rhs(1:nfr) = - g(ifr(1:nfr))
      rhs(nfr+1:lss%matrix%n) = 0.0d0
      
      call lsssol(lss,lss%matrix%n,rhs)

      d(1:nfr) = rhs(1:nfr)

      write(*,*) 'dsupn = ',maxval( abs( d(1:nfr) ) ),' deucn = ',norm2( d(1:nfr) )

      if ( norm2( d(1:nfr) ) .gt. eta * max( 1.0d0, norm2( x(1:nfr) ) ) ) then
         write(*,*) 'The norm of the direction is too big; so we skip it and increase sigma.'
         if ( 10.0d0 * sgmtrial .le. lrgstreal ) then
            sgmtrial = 10.0d0 * sgmtrial
            go to 210
         else
            write(*,*) 'sigma is already too big.'
         end if
      end if

300   continue
      
      ! Compute maximum feasible step

      amaxv(1:nfr) = lrgstreal
      
      where ( d(1:nfr) .lt. 0.0d0 .and. lind(ifr(1:nfr)) )
         amaxt(1:nfr) = 'L'
         amaxv(1:nfr) = ( lbnd(ifr(1:nfr)) - x(ifr(1:nfr)) ) / d(1:nfr)
      end where

      where ( d(1:nfr) .gt. 0.0d0 .and. uind(ifr(1:nfr)) )
         amaxt(1:nfr) = 'U'
         amaxv(1:nfr) = ( ubnd(ifr(1:nfr)) - x(ifr(1:nfr)) ) / d(1:nfr)
      end where

      amax = minval( amaxv(1:nfr) )

      write(*,*) 'amax = ',amax

      if ( amax .lt. lrgstreal ) then
         rbdnnz = count( amaxv(1:nfr) .eq. amax )
         rbdind(1:rbdnnz) = pack( (/ (i,i=1,nfr) /), amaxv(1:nfr) .eq. amax )
         rbdtype(1:rbdnnz) = pack( amaxt(1:nfr), amaxv(1:nfr) .eq. amax )
         
         write(*,*) 'rbdnnz = ',rbdnnz,' rbdind = ',rbdind(1:rbdnnz),' rbdtype = ',rbdtype(1:rbdnnz)
      else
         rbdnnz = 0
         
         write(*,*) 'rbdnnz = ',rbdnnz
      end if
      
      
      ! Test projected unitary step
      
      if ( amax .lt. 1.0d0 ) then
         write(*,*) 'Unitary step is infeasible. Trying projected unitary step.'
         
         xtrial(1:n) = x(1:n)
         xtrial(ifr(1:nfr)) = xtrial(ifr(1:nfr)) + d(1:nfr)

         call projinface(n,lind,lbnd,uind,ubnd,xtrial,nfr,ifr)
        
         call evalf(n,xtrial,ftrial,ierr,pdata)
         if ( ierr .ne. 0 ) return
    
         write(*,*) 'ftrial = ',ftrial

         if ( ftrial .lt. fbtrial ) then
            fbtrial = ftrial
            xbtrial(1:n) = xtrial(1:n)
         end if
         
         if ( ftrial .le. ftarget ) then
            write(*,*) 'ftrial <= ftarget'
            return

         else if ( ftrial .le. f ) then
            write(*,*) 'Simple decrease was obtained and since we are in the frontier, this is enough.'
            extrecom   = .true.
            infrontier = .true.
            return

         else if ( ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. sgmtrial .eq. 0.0d0 .and. &
              norm2( xtrial(ifr(1:nfr)) - x(ifr(1:nfr)) ) .le. macheps14 * norm2( x(ifr(1:nfr)) ) ) then
            if ( navoidcy .lt. 2 ) then
               navoidcy = navoidcy + 1
               write(*,*) 'The unitary step is infeasible and its projection did NOT give simple decrease.'
               write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
               write(*,*) 'Therefore, we accept it anyway.'
               return
            end if
         end if
      end if

      ! Backtracking

      gtd = dot_product( g(ifr(1:nfr)), d(1:nfr) )

      write(*,*) 'gtd = ',gtd
      
      t = min( 1.0d0, amax )

310   continue
      
      xtrial(1:n) = x(1:n)
      xtrial(ifr(1:nfr)) = xtrial(ifr(1:nfr)) + t * d(1:nfr)

      if ( t .eq. amax ) then
         where ( rbdtype(1:rbdnnz) .eq. 'L' )
            xtrial(ifr(rbdind(1:rbdnnz))) = lbnd(ifr(rbdind(1:rbdnnz)))
         elsewhere
            xtrial(ifr(rbdind(1:rbdnnz))) = ubnd(ifr(rbdind(1:rbdnnz)))
         end where
      end if
      
      call evalf(n,xtrial,ftrial,ierr,pdata)
      if ( ierr .ne. 0 ) return
    
      write(*,*) 'in backtracking, t = ',t,' ftrial = ',ftrial

      if ( ftrial .lt. fbtrial ) then
         fbtrial = ftrial
         xbtrial(1:n) = xtrial(1:n)
      end if
      
       if ( t .eq. 0.0d0 ) then
          write(*,*) 'null step in backtracking'
          xtrial(1:n) = x(1:n)
          ftrial = f
          return
          
      else if ( ftrial .le. ftarget ) then
         write(*,*) 'ftrial <= ftarget'
         return

      else if ( t .eq. amax .and. ftrial .le. f ) then
         write(*,*) 'Simple decrease was obtained and since we are in the frontier, this is enough.'
         extrecom   = .true.
         infrontier = .true.
         return

      else if ( ftrial .le. f + t * gamma * gtd ) then
         write(*,*) 'Armijo holds'
         if ( t .eq. min( 1.0d0, amax ) ) then
            extrecom   = .true.
            infrontier = .false.
         end if
         return

      else if ( t .eq. min( 1.0d0, amax ) .and. sgmtrial .eq. 0.0d0 .and. &
           ftrial .le. f + macheps14 * max( 1.0d0, abs( f ) ) .and. &
           norm2( xtrial(ifr(1:nfr)) - x(ifr(1:nfr)) ) .le. macheps14 * norm2( x(ifr(1:nfr)) ) ) then
         if ( navoidcy .lt. 2 ) then
            navoidcy = navoidcy + 1
            write(*,*) 'The unitary or the maximum step did NOT give sufficient decrease.'
            write(*,*) 'However, the Euclidian norm of the step is small and the change in f is negligible.'
            write(*,*) 'Therefore, we accept it anyway.'
            return
         end if
      end if

      ttmp = ( - gtd * t ** 2 ) / ( 2.0d0 * ( ftrial - f - t * gtd ) )
      
      if ( .not. ( tau1 * t .le. ttmp .and. ttmp .le. tau2 * t ) ) then
         t = 0.5d0 * t
      else
         t = ttmp
      end if

      go to 310

    end subroutine newtonls

  end subroutine gencanb

  
  ! *****************************************************************
  ! *****************************************************************

  subroutine projinface(n,lind,lbnd,uind,ubnd,x,nind,ind)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n,nind

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    integer, intent(in) :: ind(nind)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)

    where ( lind(ind(1:nind)) )
       x(ind(1:nind)) = max( lbnd(ind(1:nind)), x(ind(1:nind)) )
    end where
  
    where ( uind(ind(1:nind)) )
       x(ind(1:nind)) = min( ubnd(ind(1:nind)), x(ind(1:nind)) )
    end where
  
  end subroutine projinface

  ! *****************************************************************
  ! *****************************************************************

  subroutine project(n,lind,lbnd,uind,ubnd,x)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: x(n)

    where ( lind(1:n) )
       x(1:n) = max( lbnd(1:n), x(1:n) )
    end where
  
    where ( uind(1:n) )
       x(1:n) = min( ubnd(1:n), x(1:n) )
    end where

  end subroutine project

end module bmgencan
