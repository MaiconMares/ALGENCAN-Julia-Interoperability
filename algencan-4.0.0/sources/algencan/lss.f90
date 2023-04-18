module lss

  use hsl_ma57_double
  
  implicit none

  type lss_type
     integer :: nnzmax
     type(zd11_type) :: matrix
     type(ma57_control) :: control
     type(ma57_factors) :: factors
  end type lss_type
  
  public :: lssini, lssend, lssfac, lsssol
  
contains
  
  ! *****************************************************************
  ! *****************************************************************

  subroutine lssini(nnzmax,lss)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: nnzmax
    type(lss_type), intent(inout) :: lss

    ! LOCAL SCALARS
    integer :: allocerr

    lss%nnzmax = nnzmax
    
    allocate(lss%matrix%row(nnzmax),lss%matrix%col(nnzmax),lss%matrix%val(nnzmax),stat=allocerr)
    
    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In lssini, allocation error.'
       stop
    end if

    call ma57_initialize(lss%factors,lss%control)
    
  end subroutine lssini

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssend(lss)

    implicit none

    ! SCALAR ARGUMENTS
    type(lss_type), intent(inout) :: lss

    ! LOCAL SCALARS
    integer :: allocerr,info
    
    call ma57_finalize(lss%factors,lss%control,info)
    
    if ( info .ne. 0 ) then
       write(*,*) 'ERROR: In lssend, ma57_finalize returned a negative info = ',info
       stop
    end if

    deallocate(lss%matrix%row,lss%matrix%col,lss%matrix%val,stat=allocerr)

    if ( allocerr .ne. 0 ) then
       write(*,*) 'ERROR: In lssend, allocation error.'
       stop
    end if

  end subroutine lssend

  ! *****************************************************************
  ! *****************************************************************

  subroutine lssfac(lss,analyse,nneigv,rankdef)

    implicit none

    ! SCALAR ARGUMENTS
    logical, intent(in) :: analyse
    logical, intent(out) :: rankdef
    integer, intent(out) :: nneigv
    type(lss_type), intent(inout) :: lss
    
    ! LOCAL SCALARS
    type(ma57_ainfo) :: ainfo
    type(ma57_finfo) :: finfo

    if ( analyse ) then
       call  ma57_analyse(lss%matrix,lss%factors,lss%control,ainfo)
          
       if ( ainfo%flag .lt. 0 ) then
          write(*,*) 'ERROR: In lssfac, ma57_analyse returned a negative ainfo%flag = ',ainfo%flag
          stop
       end if
          
       if ( ainfo%flag .eq. 1 .or. ainfo%flag .eq. 3 ) then
          write(*,*) 'ERROR: In lssfac, ma57_analyse said there are indices out of range, ainfo%flag = ',ainfo%flag
          stop
       end if
    end if
    
    call  ma57_factorize(lss%matrix,lss%factors,lss%control,finfo)

    nneigv = finfo%neig

    rankdef = finfo%flag .eq. 4 .or. finfo%flag .eq. -5

  end subroutine lssfac

  ! *****************************************************************
  ! *****************************************************************

  subroutine lsssol(lss,n,b)

    implicit none

    ! SCALAR ARGUMENTS
    integer, intent(in) :: n
    type(lss_type), intent(inout) :: lss

    ! ARRAY ARGUMENTS
    real(kind=8), intent(inout) :: b(n)

    ! LOCAL SCALARS
    type(ma57_sinfo) :: sinfo

    call ma57_solve(lss%matrix,lss%factors,b,lss%control,sinfo)

    if ( sinfo%flag .ne. 0 ) then
       write(*,*) 'ERROR: In lssfac, ma57_solve returned a non-null sinfo%flag = ',sinfo%flag
       stop
    end if
          
  end subroutine lsssol

end module lss

