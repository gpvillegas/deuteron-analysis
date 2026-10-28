module response_functions_mod
  ! ---------------------------------------------------------------------- 
  ! all structure functions are stored in common /structure/
  !
  !             i_ang,  i_enp,  i_qcm, numdat_l
  !               i_x3,   i_x1,    i_x2, numdat_l  
  
  !   max_x1 and max_x2 are currently limited to 100 in the routine splie2.f
  !   if larger values are needed the parameter nnn needs to be increased accordingly
  
  !  max_x3 currently limited to 1000 in the routiine spline.f
  !  if larger values are needed the parameter nmax needs to be increased accordingly
  
  integer,parameter :: max_x3 = 100    ! maximal number of x3 values
  integer,parameter :: max_x1 = 75    ! maximal number of x1 values
  integer,parameter :: max_x2 = 75    ! maximal numer of x2 values


  ! tolerances fr=or identical x3, x1, x2 values
  real, parameter :: eps_x3 = 1e-5   
  real, parameter :: eps_x1 = 1e-5
  real, parameter :: eps_x2 = 1e-3
  
  real,dimension(max_x3,max_x1,max_x2) :: Rl, Rt, Rlt, Rtt 
  real,dimension(max_x3,max_x1,max_x2) :: x3
  
  ! double precision version
  
  real(kind = 8),dimension(max_x3,max_x1,max_x2):: d_Rl, d_Rt, d_Rlt, d_Rtt
  
  real(kind = 8),dimension(max_x3,max_x1,max_x2)::d_th, d_x3
  
  real,dimension(max_x1) :: x1
  real,dimension(max_x2) :: x2
  
  real(kind = 8),dimension(max_x1) :: d_x1
  real(kind = 8),dimension(max_x2) :: d_x2
  
  ! limits in x3 x3_limits(i,j,1): minimum, x3_limits(i,j,2): maximum
  real(kind = 8),dimension(max_x1, max_x2, 2) :: x3_limits 
  
  ! interpolation arrays, used internally
  
  ! calculated second derivatives with resp. to x3
  ! saved between calls
  real(kind = 8), dimension(max_x3,max_x1,max_x1), save :: d_fl2_x3, d_ft2_x3, d_flt2_x3, d_ftt2_x3
  real(kind = 8), dimension(max_x3) :: d_x3_a

  
  
  integer, dimension(max_x1,max_x2) :: numdat
  
  integer :: nset, n_x3, n_x1, n_x2
  
  ! switches for debugging and input control
  logical ::  switch_x1_x2 = .false.
  logical::  make_copy = .false.
  logical :: debug = .false.
  logical :: initialized  = .false.  !  no data has been read not ready for interpolation

  
  
  ! ---------------------------------------------------------------------- 
contains
  
  subroutine resp_read_data(iunit)
    !----------------------------------------------------------------------
    ! read structure function file and store values in common
    !
    ! these data files are usually created by the programs: short_edit, short_edit1
    ! from files in the directory h2dat/deep/
    !
    ! use free form read
    !----------------------------------------------------------------------
    
    
    implicit none
    
    real:: x1_old, x2_old, x1_l, x2_l
    
    ! variables to be interpolated
    integer:: i_x3,   i_x1,    i_x2, numdat_l
    integer:: numdat_tot
    
    integer:: j, iunit,  jset_l,i_read_status
    
    character(len = 80) :: d_name_l
    
    numdat_tot = 0
    
    if (make_copy) then
       open(66, file = 'int_data.dat', status = 'unknown')
    endif
    
    ! read all data sets
    
    !----------------------------------------------------------------------
    if (debug) write (6,'(/1x,a)') 'Contents of data read in by read_data : '
    !----------------------------------------------------------------------
    
    x1_old = -1.
    x2_old = -1.
    
    i_x3 = 0
    i_x1 = 0
    i_x2 = 0
    j = 0
    
    do
       j = j + 1
       read (iunit, * , iostat = i_read_status) numdat_l ! numer of data per set
       
       if(IS_IOSTAT_END(i_read_status)) exit   ! exit loop when EOF is encountered
       
       
       ! per set
       !----------------------------------------------------------------------
       if (debug) write (6,'(1x,a,i5)') 'read data set : ', j
       !----------------------------------------------------------------------
       if (make_copy) write(66,*) numdat_l
       !----------------------------------------------------------------------
       !
       !  if x2 and x1 were interchanged they can be switched here
       !     
       if (switch_x1_x2) then
          read (iunit, *) x2_l, x1_l
       else
          read (iunit, *) x1_l, x2_l
       endif
       
       !----------------------------------------------------------------------
       if (make_copy) write(66,*) x1_l, x2_l
       !----------------------------------------------------------------------
       
       !----------------------------------------------------------------------------
       ! update indices: this works only correcly for sorted input sequence !
       !                 => use shell scripts to sort input files
       !
       !     j: data set number
       !     d_name : data set name (character string)         
       ! data file stucture
       !  numdat_l
       !   x1(1) x2(1)
       !   j, d_name         
       !   x3(1..numdat_l).. R data ..
       !   numdat_l         
       !   x1(2) x2(1)
       !   j, d_name  
       !   x3(1...numdat_l).. R data ..
       !   numdat_l
       !   x1(3) x2(1)
       !   j, d_name
       !   x3(1...numdat_l).. R data ..
       !     ..
       !   numdat_l
       !   x1(1) x2(2)
       !   j, d_name  
       !   x3(1...numdat_l).. R data ..
       !   numdat_l
       !   x1(2) x2(2)
       !   j, d_name  
       !   x3(1...numdat_l).. R data ..
       !  ...
       !  !!! the values of x1 and x2 have to be sorted !!!
       !----------------------------------------------------------------------------
       
       if (abs(x1_l - x1_old) > eps_x1) then   ! new x1 value
          i_x1 = i_x1 + 1
          x1_old = x1_l
          x1(i_x1)    = x1_l
          d_x1(i_x1) = dble(x1_l)
       endif
       if (abs(x2_l - x2_old) > eps_x2) then    ! new x2 value
          i_x2 = i_x2 + 1
          i_x1 = 1                    ! reset x1 counter
          x2_old = x2_l
          x2(i_x2)    = x2_l
          d_x2(i_x2) = dble(x2_l)
       endif
       
       !----------------------------------------------------------------------
       if (debug) then
          write (6,'(1x,a,i5,1x, i5)') ' i_x1, i_x2 : ', i_x1, i_x2
          write (6,'(1x,a,2(1pe10.3,1x))') ' x1, x2 : ', x1_l, x2_l
       endif
       !----------------------------------------------------------------------
       numdat_tot = numdat_tot + numdat_l
       numdat(i_x1, i_x2) = numdat_l
       !
       read (iunit, *)  jset_l, d_name_l  ! d_name data set name and number, jset is ignored
       
       !----------------------------------------------------------------------
       if (debug) then
          write (6,'(a, i5,5x, a)')  'Data set # : ', j, &
               'Calculation type : '//d_name_l, & ! data set name
               'numdat_l, numdat_tot : ', numdat_l, numdat_tot
          write (47, *) numdat_l, numdat_tot
       endif
       !----------------------------------------------------------------------
       !----------------------------------------------------------------------
       if (make_copy) write(66,*) jset_l, d_name_l 
       !----------------------------------------------------------------------
       
       nset = j               ! count the number of data sets in file
       ! d_name(j) = d_name_l   ! save current data set name
       
       do i_x3 = 1, numdat_l  
          !
          ! note this format depends on what machine the file has been greated
          !
          read (iunit, *) &
               x3(i_x3, i_x1, i_x2), &       ! x3
               Rl(i_x3, i_x1, i_x2), &       ! Rl 
               Rt(i_x3, i_x1, i_x2), &       ! Rt 
               Rlt(i_x3, i_x1, i_x2), &      ! Rlt
               Rtt(i_x3, i_x1, i_x2)         ! Rtt
          
          !----------------------------------------------------------------------
          if (make_copy) &
               write(66,*) &
               x3(i_x3, i_x1, i_x2), &      ! x3
               Rl(i_x3, i_x1, i_x2), &       ! Rl 
               Rt(i_x3, i_x1, i_x2), &       ! Rt 
               Rlt(i_x3, i_x1, i_x2), &      ! Rlt
               Rtt(i_x3, i_x1, i_x2)        ! Rtt
          !----------------------------------------------------------------------
          
          
          
          !----------------------------------------------------------------------
          if (debug) then
             write (6,'(6(1pe10.3,1x))')& 
                  x3(i_x3, i_x1, i_x2), &  ! 
                  Rl(i_x3, i_x1, i_x2), &  ! Rl
                  Rt(i_x3, i_x1, i_x2), &  ! Rt
                  Rlt(i_x3, i_x1, i_x2), & ! Rlt
                  Rtt(i_x3, i_x1, i_x2)   ! Rtt
          endif
          !----------------------------------------------------------------------
          
       enddo
    enddo
    write (6,'(1x,a,i5)') ' number of data sets read     : ', nset
    
    ! fill double precision arrays
    ! create the double precision arrays
    d_x3  = dble(x3)       
    d_Rl =  dble(Rl) 
    d_Rt =  dble(Rt) 
    d_Rlt = dble(Rlt) 
    d_Rtt=  dble(Rtt)
    
    n_x3 = numdat_l
    n_x2 = i_x2
    n_x1 = i_x1
    
    write (6,'(1x,a,i5,a, i5)') ' number of x3 values  between : ',&
         minval(numdat(1:n_x1, 1:n_x2)), ' and ', maxval(numdat(1:n_x1, 1:n_x2))
    write (6,'(1x,a,i5)') ' number of x2 Values  : ', n_x2
    write (6,'(1x,a,i5)') ' number of x1 values  : ', n_x1
    write (6,'(1x,a,i10)') ' total number of grid points  : ',  numdat_tot
    
    
    
    !----------------------------------------------------------------------
    if (make_copy) then
       print *, '---------------------- copy --------'
       print*, ' close copy! '
       close(66)
    endif
    ! close input file
    close(iunit)
    
    return
  end subroutine resp_read_data
  
  function resp_initialize(file_name) result(success)
    implicit none
    ! first read all the data and print set numbers and calculation type
    
    !----------------------------------------------------------------------
    !  INITIALIZE
    !----------------------------------------------------------------------
    
    
    ! calculate derivatives for the x3 and store them in arrays. 
    ! In future calls these derivatives dont have to be calculated anymore
    
    !----------------------------------------------------------------------
    ! calculate the 1st derivative for the 1st and last point with respective
    ! to the x3
    
    ! use a simple linear estimate for derivatives at edges
    !----------------------------------------------------------------------
    ! SPLINE
    !----------------------------------------------------------------------
    character(len = *), intent(in) :: file_name  ! response function file name (full path)
    integer :: success                         ! flag for success


    
    !  local variables

    integer :: j,k,m, iq, ie, i
    integer :: i_x1, i_x2, n_x3_l


    ! file parameters
    integer, parameter :: iunit = 18
    character(len=132) :: msg
    integer:: status
    
    ! variables for spline interpolation
    real(kind = 8) :: dfl1, dft1, dflt1, dftt1
    real(kind = 8) :: dfl2, dft2, dflt2, dftt2
    real(kind = 8) :: dx3

    
    real(kind = 8), dimension(max_x3) :: l_fl, l_fl2, l_ft, l_ft2, &
         l_flt, l_flt2, l_ftt, l_ftt2

    real(kind = 8) :: dx3_a



    ! open response function file
    
    open(unit = iunit, file = file_name, status="old", iomsg=msg, iostat=status)
    
    if (status > 0) then
       write(*,*) "error opening file:", msg
       write(*,*) "i/o status:", status
       success = status
       return
    end if

    success = 0
    
    ! read data
    call resp_read_data(iunit)

    ! setup permanent interpolarion array
    
    do i_x1=1, n_x1         ! loop over x1
       do i_x2 = 1, n_x2    ! loop over x2
          n_x3_l = numdat(i_x1, i_x2) ! get number of values in x3 (can vary)
          ! if (debug) print*, ' n_x3_l = ',  n_x3_l
          do i =1 , n_x3_l    
             d_x3_a(i) = d_x3(i, i_x1, i_x2)  ! get current x3 value array
          enddo

          do i = 1, n_x3_l  ! loop over x3 values
             l_fl(i) =  d_Rl(i,i_x1,i_x2)  ! get response functions at each data pair i_x1, i_x2 
             l_ft(i) =  d_Rt(i,i_x1,i_x2)  ! for each value of x3
             l_flt(i) = d_Rlt(i,i_x1,i_x2)
             l_ftt(i) = d_Rtt(i,i_x1,i_x2)
          enddo
          m = n_x3_l
          dx3_a = d_x3_a(2) - d_x3_a(1)
          
          ! initialize fl if not natural splines are used
          dfl1 = (l_fl(2) - l_fl(1)) /dx3_a
          dfl2 = (l_fl(m) - l_fl(m-1)) / dx3_a
          
          ! these are two possible ways to intitalize the spline:
          ! see Numerical Recipes (Fortran) Vol. 1. p. 88 routine spline.
          !
          ! the result is the same, and actually the routines do exactly the same
          ! namely set the second derivative of the function at the boundary to 0.
          !
          
          !               call spline (d_x3_a,l_fl, m, dfl1, dfl2, l_fl2)    ! calculate the 2nd derivatives at each x3 value for fl
          call spline (d_x3_a,l_fl, m, 1.d30, 1.d30, l_fl2)   ! natural spline
          
          ! initialize ft
          dft1 = (l_ft(2) - l_ft(1)) /dx3_a
          dft2 = (l_ft(m) - l_ft(m-1)) / dx3_a
          !               call spline (d_x3_a,l_ft, m, dft1, dft2, l_ft2)  
          call spline (d_x3_a,l_ft, m,  1.d30, 1.d30, l_ft2)   ! calculate the 2nd derivatives at each x3 value for ft
          
          ! initialize flt
          dflt1 = (l_flt(2) - l_flt(1)) /dx3_a
          dflt2 = (l_flt(m) - l_flt(m-1)) / dx3_a
          !               call spline (d_x3_a,l_flt, m, dflt1, dflt2, l_flt2)
          call spline (d_x3_a,l_flt, m,  1.d30, 1.d30, l_flt2)  ! calculate the 2nd derivatives at each x3 value for flt
          
          ! initialize ftt
          dftt1 = (l_ftt(2) - l_ftt(1)) /dx3_a
          dftt2 = (l_ftt(m) - l_ftt(m-1)) / dx3_a
          !               call spline (d_x3_a,l_ftt, m, dftt1, dftt2, l_ftt2)
          call spline (d_x3_a,l_ftt, m,  1.d30, 1.d30, l_ftt2) !  ! calculate the 2nd derivatives at each x3 value for ftt
          
          
          do i = 1, n_x3_l                  ! store the parameters in the 3d arrays
             d_fl2_x3(i,i_x1,i_x2)   = l_fl2(i)      
             d_ft2_x3(i,i_x1,i_x2)   = l_ft2(i)
             d_flt2_x3(i,i_x1,i_x2)  = l_flt2(i)
             d_ftt2_x3(i,i_x1,i_x2)  = l_ftt2(i)
          enddo
       enddo
    enddo
    initialized = .true.
    return
  end function resp_initialize
  
  
  subroutine resp_interp(x3_v, x1_v, x2_v, fl,ft,flt, ftt, status)  
    
    !----------------------------------------------------------------------
    ! returns the response function for a given enp, qcm and angle
    ! interpolated using a cubic spline in 3 dimensions
    !
    !     input: init  if true => read new file
    !            debug if true => make lots of printouts
    !            swicth12 if true => use of x2 changes faster than x1      
    !
    !       x3_v, x1_v, x2_v : values for which the response functions
    !                               are to be calculated
    !
    ! output:  f... : response functions
    !
    ! data file is read is init = true, after reading it init is set to false
    !
    ! this is a save version where all arrays get copied 
    !----------------------------------------------------------------------
    
    implicit none
    
    real, intent(in) ::  x3_v, x1_v, x2_v
    ! output
    real, intent(out) :: fl, ft, flt, ftt
    integer, intent(out) :: status
    

    ! local arrays for x3 interpolated values
    real(kind = 8), dimension(max_x3) :: fl_x3, ft_x3, flt_x3, ftt_x3

    ! arrays of second derivatives for x3 interpolated values
    real(kind = 8), dimension(max_x3) :: fl2_x3, ft2_x3, flt2_x3, ftt2_x3
    
    ! local arrays for x2 interpolated values
    real(kind = 8), dimension(max_x2) :: fl_x2, flt_x2, ft_x2, ftt_x2 
    
    ! arrays of second derivatives for x2 interpolated values
    real(kind = 8), dimension(max_x2) :: fl2_x2, flt2_x2, ft2_x2, ftt2_x2 

    ! local arrays for x1 interpolated values
    real(kind = 8), dimension(max_x1) :: fl_x1, flt_x1, ft_x1, ftt_x1 
    
    ! arrays of second derivatives for x1 interpolated values
    real(kind = 8), dimension(max_x1) :: fl2_x1, flt2_x1, ft2_x1, ftt2_x1 
    
    ! result of final interpolation
    real(kind = 8) :: d_fl, d_flt, d_ft, d_ftt 
    
    !  local variables
    integer :: j,k,m, iq, ie, i
    integer :: i_x1, i_x2, n_x3_l
    integer :: i_x1_ok, i_x2_ok, i_x3_ok


    real(kind = 8), dimension(max_x2) :: d_x2_loc
    real(kind = 8), dimension(max_x1) :: d_x1_loc


    real(kind = 8) :: dx3_a

    logical out_of_bounds_x1, out_of_bounds_x2, out_of_bounds_x3
    
    
    
    
    !----------------------------------------------------------------------
    ! for each x1 and x2 value pair calculate the x3 interpolated values
    !    interpolate and use that the arrays are stored column wise 
    !    use spline interpolation:
    !
    !----------------------------------------------------------------------
    !    SPLINT 
    !----------------------------------------------------------------------

    i_x2_ok = 0
    do i_x1 = 1, n_x1
       i_x3_ok = 0
       do i_x2 = 1, n_x2
          n_x3_l = numdat(i_x1, i_x2) ! get number of values in x3 (can vary)
          do i = 1, n_x3_l
             fl_x3(i)    =  d_Rl(i,i_x1,i_x2)      ! get the responses for each point in x1, x2
             fl2_x3(i)   =    d_fl2_x3(i,i_x1,i_x2)
             ft_x3(i)    =  d_Rt(i,i_x1,i_x2)  
             ft2_x3(i)   =    d_ft2_x3(i,i_x1,i_x2)  
             flt_x3(i)   =  d_Rlt(i,i_x1,i_x2)  
             flt2_x3(i)  =   d_flt2_x3(i,i_x1,i_x2)  
             ftt_x3(i)   = d_Rtt(i,i_x1,i_x2) 
             ftt2_x3(i)  =   d_ftt2_x3(i,i_x1,i_x2)  
             
          enddo
          do i =1 , n_x3_l    
             d_x3_a(i) = d_x3(i, i_x1, i_x2)  ! get current x3 value array
          enddo
          m = n_x3_l
          ! check if x3_v is withing the current range of x3 grid values
          out_of_bounds_x3 =  (dble(x3_v) < d_x3_a(1)) .or. (dble(x3_v) > d_x3_a(m)) 
          if (.not. out_of_bounds_x3) then
             ! count the good x3 values
             ! and store the corresponding x1, x2 indices
             i_x3_ok = i_x3_ok + 1
             
             call splint(d_x3_a, fl_x3, fl2_x3, m, dble(x3_v), &       ! interpolate along x3 at each x1, x2 point
                  fl_x2(i_x3_ok) )                                   ! store the result in the arrauys da_fl, da_flt etc.
             
             call splint(d_x3_a, flt_x3, flt2_x3, m, dble(x3_v), &
                  flt_x2(i_x3_ok) )
             
             call splint(d_x3_a, ft_x3, ft2_x3, m, dble(x3_v),  &
                  ft_x2(i_x3_ok) )
             
             call splint(d_x3_a, ftt_x3, ftt2_x3, m, dble(x3_v), & 
                  ftt_x2(i_x3_ok) )
             d_x2_loc(i_x3_ok)  = d_x2(i_x2)
          endif
       enddo
       if (i_x3_ok .eq. 0) then
          ! not x3 data points found in range, no point to continue
          ! cannot interpolate
          status = 3
          fl = -1.
          flt = -1.
          ft = -1.
          ftt = -1.
          return
       endif
       ! check that values are not out of bounds
       ! check if x2_v is withing the current range of x2 grid values
       out_of_bounds_x2 =  (dble(x2_v) < d_x2_loc(1)) .or. (dble(x2_v) > d_x2_loc(i_x3_ok))
       if (.not. out_of_bounds_x2) then
          i_x2_ok = i_x2_ok + 1

          ! interpolate along da_fl etc. contain the x3 interpolated value for the current values ix_2 = 1,..., i_x3 = const
          ! d_x2_loc : is the array of corresonding x2 values
          
          ! setup spline for da_fij along d_x_2_loc

          m = i_x3_ok

          
          call spline (d_x2_loc,fl_x2,  m, 1.d30, 1.d30, fl2_x2)   ! natural spline
          call spline (d_x2_loc,flt_x2, m, 1.d30, 1.d30, flt2_x2)
          call spline (d_x2_loc,ft_x2,  m, 1.d30, 1.d30, ft2_x2)
          call spline (d_x2_loc,ftt_x2, m, 1.d30, 1.d30, ftt2_x2)
          
          ! interpolate along x2
          call splint(d_x2_loc, fl_x2, fl2_x2, m, dble(x2_v), &       ! interpolate along x2 at each x1 point
               fl_x1(i_x2_ok) )                                          ! store the result in the arrauys da_fl, da_flt etc.

          call splint(d_x2_loc, flt_x2, flt2_x2, m, dble(x2_v), &       ! interpolate along x2 at each x1 point
               flt_x1(i_x2_ok) )                                          ! store the result in the arrauys da_fl, da_flt etc.
          
          call splint(d_x2_loc, ft_x2, ft2_x2, m, dble(x2_v), &        ! interpolate along x2 at each x1 point
               ft_x1(i_x2_ok) )                                           ! store the result in the arrauys da_fl, da_flt etc.
          
          call splint(d_x2_loc, ftt_x2, ftt2_x2, m, dble(x2_v), &      ! interpolate along x2 at each x1 point
              ftt_x1(i_x2_ok) )                                          ! store the result in the arrauys da_fl, da_flt etc.
          d_x1_loc(i_x2_ok)  = d_x1(i_x1)
       endif
    enddo
    if (i_x2_ok .eq. 0) then
       ! not x2 data points found in range, no point to continue
       ! cannot interpolate
       status = 2
       fl = -1.
       flt = -1.
       ft = -1.
       ftt = -1.
       return
    endif

    
    ! final interpolation along x1
    ! check boundary
    out_of_bounds_x1 =  (dble(x1_v) < d_x1_loc(1)) .or. (dble(x1_v) > d_x1_loc(i_x2_ok))
    if (out_of_bounds_x1) then
       status = 1
       fl = -1.
       flt = -1.
       ft = -1.
       ftt = -1.
    else
       ! setup spline for da_fij along d_x_2_loc
       status = 0

       m = i_x2_ok
       
       call spline (d_x1_loc, fl_x1,  m, 1.d30, 1.d30, fl2_x1)   ! natural spline
       call spline (d_x1_loc, flt_x1, m, 1.d30, 1.d30, flt2_x1)
       call spline (d_x1_loc, ft_x1,  m, 1.d30, 1.d30, ft2_x1)
       call spline (d_x1_loc, ftt_x1, m, 1.d30, 1.d30, ftt2_x1)    
       ! perform the lasy interpolation in x1

       ! interpolate along x2
       call splint(d_x1_loc, fl_x1, fl2_x1, m, dble(x1_v), &       ! interpolate along x2 at each x1 point
            d_fl)                                          ! store the result in the arrauys da_fl, da_flt etc.
       
       call splint(d_x1_loc, flt_x1, flt2_x1, m, dble(x1_v), &       ! interpolate along x2 at each x1 point
            d_flt)                                          ! store the result in the arrauys da_fl, da_flt etc.
       
       call splint(d_x1_loc, ft_x1, ft2_x1, m, dble(x1_v), &        ! interpolate along x2 at each x1 point
            d_ft)                                           ! store the result in the arrauys da_fl, da_flt etc.
       
       call splint(d_x1_loc, ftt_x1, ftt2_x1, m, dble(x1_v), &      ! interpolate along x2 at each x1 point
            d_ftt)                                          ! store the result in the arrauys da_fl, da_flt etc.
       
       fl = sngl(d_fl)
       flt = sngl(d_flt)
       ft = sngl(d_ft)
       ftt = sngl(d_ftt)
    endif
    
    return 
  end subroutine resp_interp

  
  function resp_calc_sigma(k_i, k_f, theta_e, p_r, phi_r, v_l, v_t, v_lt, v_tt, s_fact) result(crs_sf)
    implicit none
    ! calculate the cross section and response functions for a set of kinematics
    ! given by k_i , k_f, theta_e (electron side)
    !         p_r, phi_r        (proton side)
    ! all angles in radians
    ! all momenta in GeV/c


    real, intent(in):: k_i, k_f, theta_e, p_r, phi_r

    ! return leptonic factors and s_fact (mott cross section and recoil factor and unit conversion)
    real, intent(out):: v_l, v_t, v_lt, v_tt, s_fact
    
    ! output
    real :: crs_sf
    
    ! local variables
    ! particle masses
    !     pm: proton mass
    !     pmn: neutron mass
    !     pmp: proton mass
    !     dm  : deuteron mass
    
    real, parameter :: pm  = 0.938279         
    real, parameter :: dm  = 1.875628
    real, parameter :: pmp = pm
    real, parameter :: pmn = 0.939565  
    real, parameter :: eb  = 0.00221
    real, parameter :: alpha = 1./137.
    real, parameter  :: pi = 3.141592653589793
    real, parameter :: dtr = pi/180.
    
    ! electron kin. variables
    real::  q2, qv, q2v
    real:: th_f, e_f, e_r, thr, thr_deg, p_f
    ! response functions
    real :: w_l, w_t, w_lt, w_tt
    integer :: i_status
    
    ! other
    real:: sfc, a_un, fc1, fc, fkin, sf0, sf
    
    real:: cthr, cthp

    
    
    !     k_i beam energy
    !     k_f  : scattered electron energy
    !     theta_e electron scattering angle
    !     p_f proton final momentum
    !     p_r  : recoil momentum
    
    !     qv  3 vector magnitude
    !     th_f proton angle (wrsp to q?)  assume relative to q
    !     e_r : energy of recoil particle (neutron)
    !     e_f : proton energy

    ! check if response functions are initialized
    if (.not. initialized) then
       print*, ' Response functions not initialized !'
       return
    endif

    ! particle energies
    e_r = sqrt(pmn**2 + p_r**2)
    e_f = (k_i - k_f) + dm - e_r  ! proton energy
    
    p_f = sqrt(e_f**2 - pm**2)
    
    ! calculate kinematics
    ! electron kinematics
    q2 = 4.*k_i*k_f* sin(theta_e/2)**2
    q2v = q2 + (k_i - k_f)**2
    qv = sqrt(q2v)
    
    ! proton kinematics
    cthr = (q2v + p_r**2 - p_f**2)/(2.*qv*p_r)
    cthp = (q2v + p_f**2 - p_r**2)/(2.*qv*p_f)

    thr = acos(cthr)
    thr_deg = thr/dtr   ! recoil angle in degrees (needed for interpolation)
    
    th_f = acos(cthp)
    
    ! calc response functions
    call resp_interp(p_r, q2, thr_deg , w_l, w_t, w_lt, w_tt, i_status) 

    ! check the intepolation status
    if (i_status .ne. 0) then
       ! bad interpolation return negative cross section with the negative of the status (for denugging)
       crs_sf = -i_status
       return
    endif
    ! return negative cross sections for unphysical kinematics
    if (abs(cthr) .gt. 1) then
       crs_sf = -10.
       return
    endif
    if (abs(cthp) .gt. 1) then
       crs_sf = -20.
       return
    endif
    
    
    v_l  = q2**2/qv**4
    v_t  = q2/2.0/qv**2 + tan(theta_e/2.0)**2
    v_tt = q2/2.0/qv**2
    v_lt = q2/qv**2*sqrt(q2/qv**2 + tan(theta_e/2.0)**2)
    
    sfc  = 4.0*k_f*k_i*cos(theta_e/2.0)**2 /3.0
    
    a_un =   sfc* (v_l*w_l + v_t*w_t + v_lt*w_lt*cos(phi_r) + v_tt*w_tt*cos(2*phi_r))
    
    !**************************************
    !       Kinematic factors
    !**************************************
    fc1 = (p_f - qv*cos(th_f))/e_r
    fc = abs(p_f/e_f + fc1)
    fkin = p_f**2 / fc
    sf0 = alpha**2/q2**2 * k_f/k_i *1.0/(2.0*dm*e_f) 
    
    !         Above derivation is from page 8 of folder
    !	sf0 = alpha**2/q2**2 * k_f/k_i *1.0/(4.0*e_r*e_f) 
    
    sf  = sf0*0.389385*1000.*1000.  ! transfer GeV/c to nb
    
    !**************************************
    !       Cross section 
    !**************************************

    s_fact = sf*fkin*sfc  ! factor by which response is multiplied to get to the v_ij w_ij terms
    
    crs_sf = sf*fkin*a_un
    
  end function resp_calc_sigma


  ! spline interpolations routines
  subroutine spline(x,y,n,yp1,ypn,y2)
    
    ! input:
    ! x, y give grid points
    ! n number of values in in x and y
    ! yp1 : first derivative at x(1) > 1e30 : set second derivatve 0 (natural spline)
    ! ypn : first deribvative at x(n) > 1e30 : set second derivatve 0 (natural spline)
    
    ! output:
    ! y2 : second derivatives
    
    
    ! min. nr. points = 3
    
    
    implicit none
    
    integer, parameter :: nmax=1000
    
    integer, intent(in) :: n
    real(kind = 8), intent(in) :: yp1, ypn
    real(kind = 8), dimension(n), intent(in) :: x(n), y(n)
    
    real(kind = 8), dimension(n), intent(out) :: y2(n)
    
    real(kind = 8), dimension(nmax) :: u
    
    integer :: i, k
    real(kind = 8) :: p, qn, un, sig
    
    if (n .gt. nmax) then
       stop ' too many points for spline'
    endif
    
    if (n .eq. 1) then
       y2(1) = 0.
       return
    endif
    
    if (yp1.gt..99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) / &
            (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
    enddo
    
    if (ypn.gt..99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    
    do k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    
    return
  end subroutine spline
  
  subroutine splint(xa,ya,y2a,n,x,y)
    
    ! input
    ! xa, ya grid data
    ! n lenght of xa and ya
    ! y2 second derivatives at xa
    
    ! x : value of x at which y should be evaluated
    
    
    ! output
    ! y : interpolated value
    
    
    implicit none
    
    integer ::  n
    real(kind = 8) ::  x, y
    real(kind = 8) :: xa(n),ya(n)
    real(kind = 8) :: y2a(n)
    
    real(kind = 8) :: h, a, b
    
    integer::  klo, khi, k
    
    klo=1
    khi=n
    
    do while (khi-klo.gt.1) 
       k=(khi+klo)/2
       if(xa(k).gt.x)then
          khi=k
       else
          klo=k
       endif
    enddo
    
    h=xa(khi)-xa(klo)
    
    if (h.eq.0.) then
       print *,'----- splint x-values ----'
       print*, xa
       print*, '----- splint y-values ----'
       print*, ya
       print*, '----- splint x  ----'
       print*, x
       stop 'bad xa input.'
    endif
    
    a = (xa(khi) - x)/h
    b = (x - xa(klo))/h
    y = a*ya(klo) + b*ya(khi) + ((a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi))*(h**2)/6.
    
    return
  end subroutine splint
  
end module response_functions_mod
