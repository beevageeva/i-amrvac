module mod_usr
  use mod_twofl


  implicit none

  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr

  double precision :: Bz0,Lx,Lz,LzN, freq,Eval,rhon00,pen00
  double precision :: rhoc00,pec00
  double precision :: rhoc01,LzC

  double precision,parameter :: ionMass=25d0

contains

  !> Read this module's parameters from a file
  subroutine usr_params_read(files)
    character(len=*), intent(in) :: files(:)
    integer                      :: n

!    namelist /usr_list/ 
!
!    do n = 1, size(files)
!       open(unitpar, file=trim(files(n)), status="old")
!       read(unitpar, usr_list, end=111)
!111    close(unitpar)
!    end do

  end subroutine usr_params_read


  subroutine dump_units()
    character(len=*), parameter :: units_filename="units.dat"
    type(arrayptr), dimension(7) :: write_vars
    integer, parameter :: I_LEN =1, I_TIME=2, I_DENS=3, I_VEL=4, I_PRES=5, I_MAGFIELD=6, I_TEMP=7

    write_vars(I_LEN)%namevar = "unit_length"
    write_vars(I_TIME)%namevar = "unit_time"
    write_vars(I_DENS)%namevar = "unit_density"
    write_vars(I_VEL)%namevar = "unit_velocity"
    write_vars(I_PRES)%namevar = "unit_pressure"
    write_vars(I_MAGFIELD)%namevar = "unit_magneticfield"
    write_vars(I_TEMP)%namevar = "unit_temperature"

    allocate(write_vars(I_LEN)%p(1), write_vars(I_TIME)%p(1), write_vars(I_DENS)%p(1), write_vars(I_VEL)%p(1), write_vars(I_PRES)%p(1), write_vars(I_MAGFIELD)%p(1), write_vars(I_TEMP)%p(1))
    write_vars(I_LEN)%p(1) = unit_length
    write_vars(I_TIME)%p(1) = unit_time
    write_vars(I_DENS)%p(1) = unit_density
    write_vars(I_VEL)%p(1) = unit_velocity
    write_vars(I_PRES)%p(1) = unit_pressure
    write_vars(I_MAGFIELD)%p(1) = unit_magneticfield
    write_vars(I_TEMP)%p(1) = unit_temperature
    call  write_formatted_file(units_filename,write_vars)

  deallocate(write_vars(I_LEN)%p, write_vars(I_TIME)%p, write_vars(I_DENS)%p, write_vars(I_VEL)%p, write_vars(I_PRES)%p, write_vars(I_MAGFIELD)%p, write_vars(I_TEMP)%p)
  end subroutine dump_units



  subroutine usr_init()
    use mod_variables
    use mod_convert 
    ! default units = 1
    unit_length        = 1d5   ! m
    unit_temperature   = 1d4   ! K
    unit_numberdensity = 3d10  ! m^-3



    usr_init_one_grid => initonegrid_usr
    usr_set_equi_vars => special_set_equi_vars
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0

    usr_set_parameters  => init_params_usr
    usr_special_bc => specialbc_usr  
    usr_internal_bc => set_vars
    usr_mask_gamma_ion_rec => set_constant_gamma_rec
    usr_mask_alpha => set_alpha
    !usr_source =>  specialsource        
    call add_convert_method2(dump_vars, 3, "jx jy jz", "_aux_") 
    call set_coordinate_system("Cartesian_2.5D")

    call twofl_activate()


  end subroutine usr_init

    subroutine set_constant_gamma_rec(ixI^L,ixO^L,w,x,gamma_ion, gamma_rec)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(inout) :: gamma_ion(ixI^S),gamma_rec(ixI^S)
  
      gamma_ion=0d0
      gamma_rec=1d0*unit_time

    end subroutine set_constant_gamma_rec


    subroutine set_alpha(ixI^L,ixO^L,w,x,alpha)
      use mod_global_parameters
      integer, intent(in) :: ixI^L, ixO^L
      double precision, intent(in) :: x(ixI^S,1:ndim)
      double precision, intent(in) :: w(ixI^S,1:nw)
      double precision, intent(inout) :: alpha(ixI^S)
  
      alpha = 2e3 * unit_time/rhon00 * exp(-((x(ixO^S,2)-xprobmin2)/LzN)**2)
    end subroutine set_alpha
  
 function dump_vars(ixI^L, ixO^L, w, x, nwc) result(wnew)
   use mod_global_parameters
   use mod_thermal_emission
   integer, intent(in)             :: ixI^L,ixO^L, nwc
   double precision, intent(in)    :: w(ixO^S, 1:nw)
   double precision, intent(in)    :: x(ixO^S,1:ndim)
   double precision    :: wnew(ixO^S, 1:nwc)
   double precision :: current(ixI^S,7-2*ndir:3)
   integer :: idirmin,idir
   call get_current(w,ixI^L,ixO^L,idirmin,current)
   wnew(ixO^S,1) =  current(ixO^S,1)
   wnew(ixO^S,2) =  current(ixO^S,2)
   wnew(ixO^S,3) =  current(ixO^S,3)
   print*,'maxcurrent', minval(wnew(ixO^S,1:3)),&
        maxval(wnew(ixO^S,1:3))
 end function dump_vars


! subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
!    use mod_global_parameters
!
!    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
!    double precision, intent(in) :: qdt, qtC, qt
!    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
!    double precision, intent(inout) :: w(ixI^S,1:nw)
!
!    double precision                :: tmp(ixO^S)
!
!    
!
!    tmp(ixO^S)=Eval*&
!        exp(-( (x(ixO^S,2)-xprobmax2) **2/Lz**2 )) * &
!        exp(-(x(ixO^S,1)**2/Lx**2 )) * cos(2d0 * dpi * freq * qt)/Bz0
!    w(ixO^S,mom_c(2)) = tmp(ixO^S) 
! end subroutine specialsource


    subroutine set_vars(level,qt,ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)

      double precision :: Efield(ixI^S)

    !print*, 'SETVARS VELy', minval(w(ixO^S,mom_c(:)))
    !print*, 'SETVARS VELy', maxval(w(ixO^S,mom_c(:)))
    
   call twofl_to_primitive(ixI^L,ixO^L,w,x)  
    w(ixI^S,rho_n_)=0d0  
    w(ixI^S,mom_n(1:ndir))=0d0  
    w(ixI^S,e_n_)=0d0  
    w(ixI^S,e_c_)=0d0 

    if(iprob==2) then
      where(x(ixO^S,2)>Lz)
        Efield(ixO^S)= Eval* &
          exp(-(x(ixO^S,1)**2/Lx**2 )) !* sin(2d0 * dpi * freq * qt)
  
        w(ixO^S,mom_c(3)) = Efield(ixO^S)/Bz0 
      endwhere
    endif
 
   call twofl_to_conserved(ixI^L,ixO^L,w,x)  
  
  end subroutine set_vars


  subroutine specialbc_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in)             :: ixG^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision                :: Efield(ixO^S)

    if(iB==3) then

      w(ixO^S,mag(2)) = 0d0
    elseif (iB==4 .and. iprob==1) then
      !Boundary condition for zmax
      call twofl_to_primitive(ixG^L,ixO^L,w,x)  
      !w(ixO^S, 1:nw) = 0d0
      Efield(ixO^S)= Eval* &
          exp(-(x(ixO^S,1)**2/Lx**2 )) !* sin(2d0 * dpi * freq * qt)
      w(ixO^S,e_c_)=0d0
      w(ixO^S,e_n_)=0d0
      w(ixO^S,rho_n_)=0d0
      w(ixO^S,mom_n(1:3))=0d0
      !w(ixO^S,mag(2)) = w(ixO^S,mag(2)) + dt * Efield(ixO^S)/dxlevel(2)
  
      w(ixO^S,mom_c(3)) = Efield(ixO^S)/Bz0 
  
      call twofl_to_conserved(ixG^L,ixO^L,w,x)  
    else
      print*, "WRONG BC DIR ", iB
    endif


  end subroutine specialbc_usr


    !> Here one can add a steady (time-independent) equi vars
    subroutine special_set_equi_vars(ixI^L,ixO^L,x,w0)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w0(ixI^S,1:number_equi_vars)

      double precision   :: n01, n02
                

      w0(ixO^S,equi_pe_n0_) = pen00
      w0(ixO^S,equi_pe_c0_) = pec00
      w0(ixO^S,equi_rho_n0_) = rhon00 
      w0(ixO^S,equi_rho_c0_) = rhoc00 * ionMass


    end subroutine special_set_equi_vars


  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
      wB0(ixO^S,2)=Bz0
      wB0(ixO^S,1)=0d0
      wB0(ixO^S,3)=0d0

  end subroutine specialset_B0

  subroutine specialset_J0(ixI^L,ixO^L,x,wJ0)
  ! Here add a time-independent background current density 
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wJ0(ixI^S,7-2*ndir:ndir)
    wJ0(ixO^S,:)=zero


  end subroutine specialset_J0




  subroutine initonegrid_usr(ixI^L,ixO^L,w,x)
    ! initialize one grid
    integer, intent(in) :: ixI^L, ixO^L
    double precision, intent(in) :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: w(ixI^S,1:nw)
    complex, parameter :: ic = dcmplx(0,1)
    double complex :: wave(ixO^S), B1, Vc, Vn,omega
    double precision :: va0, k

    w(ixO^S,1:nw) = 0d0
    if(mype .eq. 0) then
        print*, "UNITMAG ", unit_magneticfield        
        print*, "UNITTIME ", unit_time      
        print*, "UNITVEL ", unit_velocity       
    endif        
    !call twofl_to_conserved(ixI^L,ixO^L,w,x)  
  end subroutine initonegrid_usr

  
  subroutine write_formatted_file(filename,vars)
    character(len=*), intent(in) :: filename
    type(arrayptr), intent(in), dimension(:) :: vars
    integer :: i,j
    integer, parameter :: reading_unit=100
    OPEN(reading_unit, file=filename)
    write(reading_unit,'(A)',advance="no") "#"
    !!this should also work if the format is 'number of vars(A)'instead of just '(A)' which puts one value at one line
    !write(reading_unit,'(A)') (vars(i)%namevar,i=1,size(vars))
    write(reading_unit,*) (vars(i)%namevar,i=1,size(vars))
    do i=1,size(vars(1)%p)
      !write(reading_unit,'(E2.5)') (vars(j)%p(i),j=1,size(vars))
      write(reading_unit,*) (vars(j)%p(i),j=1,size(vars))
    enddo
    CLOSE(reading_unit)
  end subroutine write_formatted_file


  subroutine init_params_usr()
    double precision :: unit_ele
    double precision, parameter :: echarge =1.6e-19 !Coulomb

    call dump_units()

    Lx = 1d0 !Lx=100km       
    !Lz = 5d-1   !used for setting source close to upper BC 
    Lz = 11d0   !lower height of the pert at the upper BC 
    LzN=0.4 !40km  
    unit_ele = unit_length  * unit_magneticfield / unit_time
    Eval = -62.4*1d-3/unit_ele 
    freq =  1d-2 * unit_time    

    Bz0=-5d-5/unit_magneticfield

    rhon00=1d3
    pen00=1d0

    rhoc00=1d0
    pec00=1d0

    rhoc01=2d0
    LzC=0.5 !50km, Case 2 

    twofl_etah=ionMass*mp_SI/echarge / (unit_magneticfield * unit_time)

    if(mype .eq. 0) then   
      print*, "DRIVER FREQ. = ", 1d0/freq
      print*, "V_ini ", Eval/Bz0 * unit_length/unit_time
      print*, "ETAH ", twofl_etah
    endif

  end subroutine init_params_usr





  function rargsort(a) result(b)
    ! Returns the indices that would sort an array.
    !
    ! Arguments
    ! ---------
    !
    double precision, intent(in):: a(:)   ! array of numbers
    integer :: b(size(a))         ! indices into the array 'a' that sort it
    !
    ! Example
    ! -------
    !
    ! rargsort([4.1_dp, 2.1_dp, 2.05_dp, -1.5_dp, 4.2_dp]) ! Returns [4, 3, 2, 1, 5]
    
    integer :: N                           ! number of numbers/vectors
    integer :: i,imin                      ! indices: i, i of smallest
    integer :: temp1                       ! temporary
    double precision:: temp2
    double precision:: a2(size(a))
    a2 = a
    N=size(a)
    do i = 1, N
        b(i) = i
    end do
    do i = 1, N-1
        ! find ith smallest in 'a'
        imin = minloc(a2(i:),1) + i - 1
        ! swap to position i in 'a' and 'b', if not already there
        if (imin /= i) then
            temp2 = a2(i); a2(i) = a2(imin); a2(imin) = temp2
            temp1 = b(i); b(i) = b(imin); b(imin) = temp1
        end if
    end do
  end function






end module mod_usr
