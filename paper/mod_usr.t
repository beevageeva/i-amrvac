module mod_usr
  use mod_twofl


  implicit none

  type arrayptr
    double precision,dimension(:), pointer :: p
    character(len=30) ::  namevar
  end type arrayptr

  double precision :: Bz0,Lx,Ly,Lz,freq,Eval


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





  subroutine usr_init()
    use mod_variables
    use mod_convert 
    ! default units = 1
    unit_length        = 1d5   ! m
    unit_temperature   = 1d3   ! K
    unit_numberdensity = 3d10  ! m^-3



    usr_init_one_grid => initonegrid_usr
    usr_set_equi_vars => special_set_equi_vars
    usr_set_B0              => specialset_B0
    usr_set_J0              => specialset_J0

    usr_set_parameters  => init_params_usr
    usr_special_bc => specialbc_usr  
    usr_internal_bc => set_vars
    !usr_source =>  specialsource        
    call add_convert_method2(dump_vars, 3, "jx jy jz", "_aux_") 
    call set_coordinate_system("Cartesian_3D")

    call twofl_activate()


  end subroutine usr_init
  
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


 subroutine specialsource(qdt,ixI^L,ixO^L,iw^LIM,qtC,wCT,qt,w,x)
    use mod_global_parameters

    integer, intent(in) :: ixI^L, ixO^L, iw^LIM
    double precision, intent(in) :: qdt, qtC, qt
    double precision, intent(in) :: x(ixI^S,1:ndim), wCT(ixI^S,1:nw)
    double precision, intent(inout) :: w(ixI^S,1:nw)

    double precision                :: tmp(ixO^S)

    w(ixO^S, 1:nw) = 0d0
    tmp(ixO^S)=-dt * 2d0 * Eval*&
        x(ixO^S,2)/Ly**2 * exp(-(x(ixO^S,2)**2/Ly**2 ))*&
        exp(-( (x(ixO^S,3)-xprobmax3) **2/Lz**2 )) * &
        sin(2d0 * dpi * freq * qt)
    !print*, 'MINMAX', minval(tmp(ixO^S)),&
   !maxval(tmp(ixO^S))   
    w(ixO^S,mag(3))=w(ixO^S,mag(3)) + tmp(ixO^S) 
    tmp(ixO^S)=Eval*exp(-(x(ixO^S,2)**2/Ly**2 ))*&
        exp(-( (x(ixO^S,3)-xprobmax3) **2/Lz**2 )) * &
        exp(-(x(ixO^S,1)**2/Lx**2 )) * cos(2d0 * dpi * freq * qt)/Bz0
    w(ixO^S,mom_c(2)) = tmp(ixO^S) 
 end subroutine specialsource


    subroutine set_vars(level,qt,ixI^L,ixO^L,w,x)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L,level
      double precision, intent(in)    :: qt
      double precision, intent(inout) :: w(ixI^S,1:nw)
      double precision, intent(in)    :: x(ixI^S,1:ndim)


    !print*, 'SETVARS VELy', minval(w(ixO^S,mom_c(:)))
    !print*, 'SETVARS VELy', maxval(w(ixO^S,mom_c(:)))
    
    !w(ixI^S,rho_n_)=0d0  
    !w(ixI^S,mom_n(1:ndir))=0d0  
    w(ixI^S,e_n_)=0d0  
    w(ixI^S,e_c_)=0d0  
   !call twofl_to_conserved(ixI^L,ixO^L,w,x)  
  
  end subroutine set_vars


  subroutine specialbc_usr(qt,ixG^L,ixO^L,iB,w,x)
    integer, intent(in)             :: ixG^L, ixO^L, iB
    double precision, intent(in)    :: qt, x(ixG^S,1:ndim)
    double precision, intent(inout) :: w(ixG^S,1:nw)

    double precision                :: Efield(ixO^S)
    double precision                :: tmp(ixO^S)



    !Boundary condition for Bz
    w(ixO^S, 1:nw) = 0d0
    Efield(ixO^S)= Eval* exp(-(x(ixO^S,2)**2/Ly**2 ))*&
        exp(-(x(ixO^S,1)**2/Lx**2 )) * sin(2d0 * dpi * freq * qt)


    tmp(ixO^S)= -2d0 * Efield(ixO^S)*&
        x(ixO^S,2)/Ly**2 
    !print*, 'MINMAX', minval(tmp(ixO^S)),&
   !maxval(tmp(ixO^S))   
    w(ixO^S,mag(3))=w(ixO^S,mag(3)) + dt * tmp(ixO^S) 

   !!!!!!!!
   ! Update horizontal components of B
   ! Spencer mod 2023/09/27

   !TODO
   !Update Bx in ghost cell
   !Bx(t+1) = Bx(t) + dt * partial Ey / partial z   
   !How to code this?
   !w(ixO^S,mag(1)) = w(ixO^S,mag(1)) + dt * !partial Ey / partial z

   !TODO
   !Update By in ghost cell
   !By(t+1) = By(t) - dt * partial Ex / partial z   
   !How to code this?
   w(ixO^S,mag(2)) = w(ixO^S,mag(2)) - dt * Efield(ixO^S)/dxlevel(3)


    tmp(ixO^S)=Efield(ixO^S)/Bz0
    w(ixO^S,mom_c(2)) = tmp(ixO^S) 

  end subroutine specialbc_usr


    !> Here one can add a steady (time-independent) equi vars
    subroutine special_set_equi_vars(ixI^L,ixO^L,x,w0)
      use mod_global_parameters
      integer, intent(in)             :: ixI^L,ixO^L
      double precision, intent(in)    :: x(ixI^S,1:ndim)
      double precision, intent(inout) :: w0(ixI^S,1:number_equi_vars)

      double precision   :: n01, n02
                

      w0(ixO^S,equi_pe_n0_) = 1d0
      w0(ixO^S,equi_pe_c0_) = 1d0
      w0(ixO^S,equi_rho_n0_) = 1d0
      w0(ixO^S,equi_rho_c0_) = 1d0


    end subroutine special_set_equi_vars


  subroutine specialset_B0(ixI^L,ixO^L,x,wB0)
  ! Here add a time-independent background magnetic field
    integer, intent(in)           :: ixI^L,ixO^L
    double precision, intent(in)  :: x(ixI^S,1:ndim)
    double precision, intent(inout) :: wB0(ixI^S,1:ndir)
      wB0(ixO^S,3)=Bz0
      wB0(ixO^S,1:2)=0d0

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
        print*, "ALPHA ", 2e3*unit_time*unit_density/(unit_numberdensity)     
    endif        
    !call twofl_to_conserved(ixI^L,ixO^L,w,x)  
  end subroutine initonegrid_usr

  


  subroutine init_params_usr()
    double precision :: unit_ele


    Lx = 1d0        
    Ly = 1d0        
    Lz = 5d-1        
    unit_ele = unit_length  * unit_magneticfield / unit_time
    Eval = 62.4*1d-3/unit_ele 
    freq =  1d-2 * unit_time       
    print*, 1d0/freq

    !Eval = Eval*1d8
    Bz0=-5d-5/unit_magneticfield
    !Bz0=Bz0*1d-3
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
