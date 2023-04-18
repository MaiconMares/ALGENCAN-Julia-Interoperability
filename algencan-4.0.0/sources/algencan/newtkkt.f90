module bmnewtkkt

  use lss, only: lss_type,lssini,lssend,lssfac,lsssol
  use memev, only: memev_type,memevini,memevend,memevscale,memevalc,memevalg,memevalj
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

  private
  
  ! SCALAR PARAMETERS
  integer, parameter :: lrgstint = huge( 1 )
  real(kind=8), parameter :: lrgstreal = huge( 1.0d0 ), macheps = epsilon( 1.0d0 ), &
       macheps12 = sqrt( macheps ), machepsinv = 1.0d0 / macheps, sgmsmall = macheps12, &
       macheps12inv = 1.0d0 / macheps12, sgmmin = macheps, sgmmax = machepsinv, &
       hmin = macheps12, hmax = macheps12inv
  
  public :: newtkkt,newtkktb

contains

  ! *****************************************************************
  ! *****************************************************************

  subroutine newtkkt(evalf,evalg,evalc,evalj,evalhl,jnnzmax,hlnnzmax, &
       n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl,epsopt, &
       maxit,scale,corrin,csupn,ssupn,nlpsupn,bdsvio,feasfound,feasx, &
       feaslambda,feascsupn,feasssupn,feasnlpsupn,feasbdsvio,infeasfound, &
       infeasx,infeaslambda,infeascsupn,infeasssupn,infeasnlpsupn,infeasbdsvio, &
       iter,ierr,istop,pdata)

    implicit none
    
    ! PROCEDURES
    procedure(fsub) :: evalf
    procedure(gsub) :: evalg
    procedure(csub) :: evalc
    procedure(jsub) :: evalj
    procedure(hlsub) :: evalhl

    ! SCALAR ARGUMENTS
    logical, intent(in) :: corrin,scale
    logical, intent(out) :: feasfound,infeasfound
    integer, intent(in) :: hlnnzmax,jnnzmax,m,maxit,n,p
    integer, intent(out) :: ierr,istop,iter
    real(kind=8), intent(in) :: epsfeas,epscompl,epsopt
    real(kind=8), intent(out) :: bdsvio,csupn,nlpsupn,ssupn,feasbdsvio,feascsupn, &
         feasnlpsupn,feasssupn,infeasbdsvio,infeascsupn,infeasnlpsupn,infeasssupn
    type(c_ptr), optional, intent(in) :: pdata

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: lambda(m+p),x(n)
    real(kind=8), intent(out) :: feaslambda(m+p),feasx(n),infeaslambda(m+p),infeasx(n)

    ! LOCAL SCALARS
    type(lss_type) :: lss
    type(memev_type) :: memev

    ierr = 0

    call lssini(hlnnzmax,lss)

    call memevini(evalf,evalg,evalc,evalj,evalhl,jnnzmax,n,m,p,memev,pdata)    

    if ( scale ) then
       call memevscale(memev,n,x,m,p,ierr)
       if ( ierr .ne. 0 ) return
    end if

    call newtkktb(lss,memev,n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas,epscompl, &
         epsopt,maxit,corrin,csupn,ssupn,nlpsupn,bdsvio,feasfound,feasx,feaslambda, &
         feascsupn,feasssupn,feasnlpsupn,feasbdsvio,infeasfound,infeasx,infeaslambda, &
         infeascsupn,infeasssupn,infeasnlpsupn,infeasbdsvio,iter,ierr,istop)
    
    if ( ierr .ne. 0 ) return

    call memevend(memev)

    call lssend(lss)
    
  end subroutine newtkkt

  ! *****************************************************************
  ! *****************************************************************

  subroutine newtkktb(lss,memev,n,x,lind,lbnd,uind,ubnd,m,p,lambda,epsfeas, &
       epscompl,epsopt,maxit,corrin,csupn,ssupn,nlpsupn,bdsvio,feasfound,feasx, &
       feaslambda,feascsupn,feasssupn,feasnlpsupn,feasbdsvio,infeasfound,infeasx, &
       infeaslambda,infeascsupn,infeasssupn,infeasnlpsupn,infeasbdsvio,iter,ierr,istop)

    implicit none
    
    ! SCALAR ARGUMENTS
    logical, intent(in) :: corrin
    logical, intent(out) :: feasfound,infeasfound
    integer, intent(in) :: m,maxit,n,p
    integer, intent(inout) :: ierr
    integer, intent(out) :: istop,iter
    real(kind=8), intent(in) :: epsfeas,epscompl,epsopt
    real(kind=8), intent(out) :: bdsvio,csupn,nlpsupn,ssupn,feasbdsvio,feascsupn, &
         feasnlpsupn,feasssupn,infeasbdsvio,infeascsupn,infeasnlpsupn,infeasssupn
    type(lss_type), intent(inout) :: lss
    type(memev_type), intent(inout) :: memev

    ! ARRAY ARGUMENTS
    logical, intent(in) :: lind(n),uind(n)
    real(kind=8), intent(in) :: lbnd(n),ubnd(n)
    real(kind=8), intent(inout) :: lambda(m+p),x(n)
    real(kind=8), intent(out) :: feaslambda(m+p),feasx(n),infeaslambda(m+p),infeasx(n)

    ! LOCAL SCALARS
    logical :: rankdef
    integer :: a,b,hlnnz,i,j,k,nelem,nfr,nkeepel,nneigv,pinc
    real(kind=8) :: sgmtrial
    
    ! LOCAL ARRAYS
    logical :: frvar(n),inc(m+p),lequ(n),workl(m+p)
    integer :: frind(n),idiag(n),incind(p),renumbered(n)
    real(kind=8) :: diagv(n),lplus(m+p),nl(n),origv(n),rhs(n+m+p),v(n),w(n),workr(m+p),xplus(n)

    iter = 0

    feasfound = .false.
    infeasfound = .false.
    
    sgmtrial = 0.0d0
    
100 continue
    
    ! Compute constraints

    call memevalc(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    if ( any( isnan( memev%c(1:m+p) ) ) ) then
       write(*,*) 'newtkktb: there is a NaN in the constraints'
       ierr = -1
       return
    end if
    
    ! Compute gradient of the Lagrangian
    
    call memevalg(n,x,ierr,memev)
    if ( ierr .ne. 0 ) return

    if ( any( isnan( memev%scaledg(1:n) ) ) ) then
       write(*,*) 'newtkktb: there is a NaN in the gradient of the objective function'
       ierr = -1
       return
    end if
    
    nl(1:n) = memev%scaledg(1:n)

    ! Selected inequality constraints' gradients are the ones that
    ! contribute to the gradient of the Lagrangian because they have a
    ! non-null multiplier plus the ones that will contribute to the
    ! KKT system because max( g(x), - lambda ) = g(x).

    if ( m + p .gt. 0 ) then
       workl(1:m) = .true.
       workl(m+1:m+p) = lambda(m+1:m+p) .ne. 0.0d0 .or. memev%scaledc(m+1:m+p) .ge. - lambda(m+1:m+p)

       if ( any( workl(1:m+p) ) ) then
          call memevalj(n,x,m,p,workl,ierr,memev)
          if ( ierr .ne. 0 ) return
          
          do j = 1,m + p
             if ( lambda(j) .ne. 0.0d0 ) then
                a = memev%jsta(j)
                b = a + memev%jlen(j) - 1

                if ( any( isnan( memev%scaledjval(a:b) ) ) .or. any( isnan( nl(1:n) ) ) ) then
                   write(*,*) 'newtkktb: there is a NaN in the Jacobian of the constraints'
                   ierr = -1
                   return
                end if

                nl(memev%jvar(a:b)) = nl(memev%jvar(a:b)) + lambda(j) * memev%scaledjval(a:b)
             end if
          end do
       end if
    end if

    ! Compute the value of the multipliers associated with bound
    ! constraints (and their effect in the gradient of the Lagrangian)

    lequ(1:n) = lind(1:n) .and. uind(1:n) .and. lbnd(1:n) .eq. ubnd(1:n)

    where( lequ(1:n) .and. nl(1:n) .le. 0.0d0 )
       v(1:n)  = 0.0d0
       w(1:n)  = - nl(1:n)
       nl(1:n) = 0.0d0
    end where

    where( lequ(1:n) .and. nl(1:n) .gt. 0.0d0 )
       v(1:n)  = nl(1:n)
       w(1:n)  = 0.0d0
       nl(1:n) = 0.0d0
    end where

    where( .not. lequ(1:n) .and. &
         ( .not. lind(1:n) .or. ( lind(1:n) .and. lbnd(1:n) .lt. x(1:n) ) ) .and. &
         ( .not. uind(1:n) .or. ( uind(1:n) .and. ubnd(1:n) .gt. x(1:n) ) ) )
       v(1:n)  = 0.0d0
       w(1:n)  = 0.0d0
    end where

    where( .not. lequ(1:n) .and. lind(1:n) .and. x(1:n) .le. lbnd(1:n) )
       v(1:n)  = nl(1:n)
       w(1:n)  = 0.0d0
       nl(1:n) = 0.0d0
    end where
    
    where( .not. lequ(1:n) .and. uind(1:n) .and. x(1:n) .ge. ubnd(1:n) )
       v(1:n)  = 0.0d0
       w(1:n)  = - nl(1:n)
       nl(1:n) = 0.0d0
    end where
    
    ! Check stopping criteria

    nlpsupn = maxval( abs( nl(1:n) ) )
    csupn   = max( 0.0d0, max( maxval( abs( memev%c(1:m) ) ), maxval( memev%c(m+1:m+p) ) ) )
    bdsvio  = max( 0.0d0, max( maxval( lbnd(1:n) - x(1:n), lind(1:n) ), maxval( x(1:n) - ubnd(1:n), uind(1:n) ) ) )
    ssupn   = max( 0.0d0, maxval( min( - memev%c(m+1:m+p), lambda(m+1:m+p) ) ), &
                          maxval( min( x(1:n) - lbnd(1:n), v(1:n) ), lind(1:n) ), &
                          maxval( min( ubnd(1:n) - x(1:n), w(1:n) ), uind(1:n) ) )
    
    write(*,*) 'newtkktb: iter = ',iter,' csupn = ',csupn,' bdsvio = ',bdsvio,' ssupn = ',ssupn,' nlpsupn = ',nlpsupn
    
    open(20,file='newtkkt-interrupted-tabline.txt')
    write(20,9000) n,m+p,csupn,ssupn,nlpsupn,bdsvio,iter
    close(20)

    ! Save the lastest feasible visited point or the best infeasible visited point
    
    if ( max( csupn, bdsvio ) .le. epsfeas ) then
       feasfound = .true.
       feasx(1:n) = x(1:n)
       feaslambda(1:m+p) = lambda(1:m+p)
       feascsupn = csupn
       feasssupn = ssupn
       feasnlpsupn = nlpsupn
       feasbdsvio = bdsvio
       
    else if ( .not. infeasfound .or. ( infeasfound .and. &
         max( csupn, bdsvio ) .le. max( infeascsupn, infeasbdsvio ) ) ) then
       infeasfound = .true.
       infeasx(1:n) = x(1:n)
       infeaslambda(1:m+p) = lambda(1:m+p)
       infeascsupn = csupn
       infeasssupn = ssupn
       infeasnlpsupn = nlpsupn
       infeasbdsvio = bdsvio
    end if
    
    if ( nlpsupn .le. epsopt .and. max( csupn, bdsvio ) .le. epsfeas .and. ssupn .le. epscompl ) then
       write(*,*) 'newtkktb: success!'
       istop = 0
       return
    end if

    if ( iter .ge. maxit ) then
       write(*,*) 'newtkktb: maximum number of iterations reached!'
       istop = 1
       return
    end if
    
    ! Iterate

    iter = iter + 1
    
    ! Set equations of the KKT system that will be considered at the
    ! current Newton step

    inc(1:m) = .true.

    ! In the case of an inequality constraint of the form g(x) <= 0
    ! with an associated Lagrange multiplier lambda >=0, we consider
    ! the equation max{ g(x), - lambda } = 0.

    ! If the maximum is attained at g(x), the constraint g(x) = 0 is
    ! considered. Otherwise, the constraint lambda = 0 is
    ! considered. In the latter case, the effect in lambda is already
    ! known (it new value will be zero) and it can be ignored.

    where( memev%scaledc(m+1:m+p) .ge. - lambda(m+1:m+p) )
       inc(m+1:m+p) = .true.
       lplus(m+1:m+p) = lambda(m+1:m+p)
    elsewhere
       inc(m+1:m+p) = .false.
       lplus(m+1:m+p) = 0.0d0
    end where

    pinc = count( inc(m+1:m+p) )

    ! We consider now the case in which g(x) is a lower (upper) bound
    ! constraint of the form li <= xi (xi <= ui).

    ! If the maximum is attained at g(x), the effect in xi is already
    ! known (it new value will be li (ui)) and it can be ignored.

    frvar(1:n) = .true.
    xplus(1:n) = x(1:n)
    
    where( lind(1:n) .and. lbnd(1:n) - x(1:n) .ge. - v(1:n) )
       frvar(1:n) = .false.
       xplus(1:n) = lbnd(1:n)
    end where
       
    where( uind(1:n) .and. x(1:n) - ubnd(1:n) .ge. - w(1:n) )
       frvar(1:n) = .false.
       xplus(1:n) = ubnd(1:n)
    end where

    nfr = count( frvar(1:n) )

    ! If nfr + m + pinc = 0 then the next iterate is already determined
    
    if ( nfr + m + pinc .eq. 0 ) then
       x(1:n) = xplus(1:n)
       lambda(m+1:m+p) = lplus(m+1:m+p)
       go to 100
    end if

    ! Set (full) jacobian

    ! Hessian of the Lagrangian
    
    workr(1:m+p) = lambda(1:m+p)

    if ( memev%scale ) then
       workr(1:m+p) = memev%sc(1:m+p) * workr(1:m+p) / memev%sf
    end if
    
    if ( memev%pdatapresent ) then
       call memev%evalhl(n,x,m,p,workr,lss%nnzmax,.true., &
            hlnnz,lss%matrix%row,lss%matrix%col,lss%matrix%val,ierr,memev%pdata)
    else
       call memev%evalhl(n,x,m,p,workr,lss%nnzmax,.true., &
            hlnnz,lss%matrix%row,lss%matrix%col,lss%matrix%val,ierr)
    end if
    
    if ( memev%scale ) then
       lss%matrix%val(1:hlnnz) = memev%sf * lss%matrix%val(1:hlnnz)
    end if

    if ( ierr .ne. 0 ) return

    if ( any( isnan( lss%matrix%val(1:hlnnz) ) ) ) then
       write(*,*) 'newtkktb: there is a NaN in the Hessian of the Lagrangian'
       ierr = -1
       return
    end if
    
    ! Remove rows and columns related to fixed variables

    renumbered(1:n) = unpack( (/ (i,i=1,nfr) /), frvar(1:n), 0 )
      
    where ( renumbered(lss%matrix%row(1:hlnnz)) .ne. 0 .and. renumbered(lss%matrix%col(1:hlnnz)) .ne. 0 )
       lss%matrix%row(1:hlnnz) = renumbered(lss%matrix%row(1:hlnnz))
       lss%matrix%col(1:hlnnz) = renumbered(lss%matrix%col(1:hlnnz))
    elsewhere
       lss%matrix%row(1:hlnnz) = 0
       lss%matrix%col(1:hlnnz) = 0
    end where

    nkeepel = count( lss%matrix%row(1:hlnnz) .ne. 0 )

    lss%matrix%ne = nkeepel
    
    lss%matrix%val(1:nkeepel) = pack( lss%matrix%val(1:hlnnz), lss%matrix%row(1:hlnnz) .ne. 0 )
      
    lss%matrix%row(1:nkeepel) = pack( lss%matrix%row(1:hlnnz), lss%matrix%row(1:hlnnz) .ne. 0 )
    lss%matrix%col(1:nkeepel) = pack( lss%matrix%col(1:hlnnz), lss%matrix%col(1:hlnnz) .ne. 0 )
      
    hlnnz = nkeepel
    
    ! Set pointers to the diagonal elements

    idiag(1:nfr) = 0
      
    do k = 1,hlnnz
       if ( lss%matrix%row(k) .eq. lss%matrix%col(k) ) then
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
          if ( hlnnz + 1 .gt. lss%nnzmax ) then
             write(*,*) 'newtkktb: lack of memory! (place 1)'
             ierr = -1
             return
          end if
          
          hlnnz = hlnnz + 1
          
          lss%matrix%row(hlnnz) = i
          lss%matrix%col(hlnnz) = i
          lss%matrix%val(hlnnz) = 0.0d0
            
          idiag(i) = hlnnz
          origv(i) = 0.0d0
          diagv(i) = 0.0d0
       end if
    end do

    ! Concatenate Jacobian of the constraints

    k = 0
    do j = 1,m + p
       if ( inc(j) ) then
          nelem = memev%jlen(j)
          
          a = memev%jsta(j)
          b = a + nelem - 1

          nkeepel = count( frvar(memev%jvar(a:b)) )
          
          if ( hlnnz + nkeepel + 1 .gt. lss%nnzmax ) then
             write(*,*) 'newtkktb: lack of memory! (place 2)'
             ierr = -1
             return
          end if

          k = k + 1
          
          lss%matrix%row(hlnnz+1:hlnnz+nkeepel) = renumbered(pack( memev%jvar(a:b), frvar(memev%jvar(a:b)) ))
          lss%matrix%col(hlnnz+1:hlnnz+nkeepel) = nfr + k
          lss%matrix%val(hlnnz+1:hlnnz+nkeepel) = pack( memev%scaledjval(a:b), frvar(memev%jvar(a:b)) )
          
          lss%matrix%row(hlnnz+nkeepel+1) = nfr + k
          lss%matrix%col(hlnnz+nkeepel+1) = nfr + k
          lss%matrix%val(hlnnz+nkeepel+1) = - sgmsmall
          
          hlnnz = hlnnz + nkeepel + 1
       end if
    end do

    if ( k .ne. m + pinc ) then
       write(*,*) 'newtkktb: there is something wrong! (k .ne. m + pinc)'
       stop
    end if
    
    lss%matrix%n = nfr + m + pinc
    lss%matrix%ne = hlnnz

200 continue
    
    if ( .not. sgmtrial .le. lrgstreal ) then
       write(*,*) 'newtkktb: sgmtrial larger than lrgstreal'
       istop = 3
       return
    end if
    
    ! Factorize

    lss%matrix%val(idiag(1:nfr)) = origv(1:nfr) + sgmtrial

    call lssfac(lss,.true.,nneigv,rankdef)

    if ( nneigv .ne. m + pinc .or. rankdef ) then
       if ( sgmtrial .eq. 0.0d0 ) then
          sgmtrial = sgmsmall * max( hmin, min( maxval( abs( diagv(1:nfr) ) ), hmax ) )
          sgmtrial = max( sgmmin, min( sgmtrial, sgmmax ) )

       else
          if ( corrin ) then
             sgmtrial = 10.0d0 * sgmtrial

          else
             write(*,*) 'newtkktb: wrong interia and not corrin'
             istop = 2
             return
          end if
          
          go to 200
       end if
    end if

    ! Set right-hand side

    frind(1:nfr) = pack( (/ (i,i=1,n) /), frvar(1:n) )
    incind(1:pinc) = pack( (/ (j,j=m+1,m+p) /), inc(m+1:m+p) )

    rhs(1:nfr) = - nl(frind(1:nfr))
    rhs(nfr+1:nfr+m) = - memev%scaledc(1:m)
    rhs(nfr+m+1:nfr+m+pinc) = - memev%scaledc(incind(1:pinc))

    ! Solve
    
    call lsssol(lss,lss%matrix%n,rhs)

    ! Update and iterate
    
    x(1:n) = xplus(1:n)
    x(frind(1:nfr)) = x(frind(1:nfr)) + rhs(1:nfr)

    lambda(1:m) = lambda(1:m) + rhs(nfr+1:nfr+m)

    lambda(m+1:m+p) = lplus(m+1:m+p)
    lambda(incind(1:pinc)) = lambda(incind(1:pinc)) + rhs(nfr+m+1:nfr+m+pinc)

    go to 100

9000 format(2(1X,I6),4(1X,1P,D7.1),6(1X,I8))

  end subroutine newtkktb

end module bmnewtkkt
