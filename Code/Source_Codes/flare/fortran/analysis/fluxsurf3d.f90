module flare_fluxsurf3d
  use iso_fortran_env
  use moose_bfit
  use flare_poincare_map
  implicit none
  private



  ! numerical parameters for flux surface construction
  integer, public :: &
     npoints = 1024, &  ! number of points used for Poincare maps
     nctrl   = 0,    &  ! number of control points for B-Spline (0: automatic refinement, see eps)
     k       = 4,    &  ! B-Spline order
     fit_method = 0     ! 0: linear, 1: non-linear

  real(real64), public :: &
     eps     = 1.d-5, & ! required accuracy for automatic refinement of B-Spline multifit
     lambda1 = 0.d0,  & ! regularization parameter for non-linear fit
     lambda2 = 1.e-5

  logical, public :: &
     knot_balancing = .true.



  ! representation of flux surface in 3D .......................................
  type, public :: fluxsurf3d
     type(poincare_map), allocatable :: section(:)
     type(bfit), allocatable :: slice(:)
     integer :: nphi

     contains
     procedure :: save, savenc
  end type fluxsurf3d


  interface fluxsurf3d
     procedure :: construct_fluxsurf3d
  end interface
  ! fluxsurf3d .................................................................



  public :: &
     loadnc_fluxsurf3d

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function construct_fluxsurf3d(r0, nsym, nphi, phi0, theta0, param, report) result(this)
  !
  ! Construct *nphi* flux surface contours with toroidal symmetry *nsym* from
  ! Poincare map for reference point *r0*.
  !
  ! Optional parameters:
  !    phi0      Lower boundary [rad] for toroidal domain [phi0, phi0 + pi2/nsym] (default: 0)
  !    theta0    Poloidal angle [rad] for lower boundary of domain.
  !    param     Select parametrization ("magnetic_angle" (default) or "arclength").
  !    report    Activate screen output.
  !
  use moose_error
  use moose_math,  only: pi2
  use moose_utils, only: progress_bar
  real(real64),     intent(in) :: r0(3)
  integer,          intent(in) :: nsym, nphi
  real(real64),     intent(in), optional :: phi0, theta0
  character(len=*), intent(in), optional :: param
  logical,          intent(in), optional :: report
  type(fluxsurf3d)             :: this

  logical :: print_parameters, print_progress_bar
  character(len=32) :: param_
  real(real64), pointer :: theta(:), rz(:,:)
  real(real64) :: phi0_, chisq_dof(-1:1)
  integer :: i, min_nctrl, max_nctrl


  print_parameters   = .false.;   if (present(report)) print_parameters   = report
  print_progress_bar = .false.;   if (present(report)) print_progress_bar = report


  ! construct flux surface geometry from Poincare map
  if (print_parameters) print 1001, npoints, r0(1:2), r0(3) / pi2 * 360.d0
  this%nphi = nphi
  phi0_ = 0.d0;   if (present(phi0)) phi0_ = phi0
  allocate (this%section(0:nphi-1), source=poincare_maps(r0, 1, phi0_, nsym, npoints, nphi))
 1001 format(3x,"- Constructing flux surface geometry from Poincare map:",/ &
             8x,"number of return points: ",i0,/ &
             8x,"reference point: (",f0.3," m, ",f0.3," m, ",f0.3," deg)"/)


  ! screen output
  if (print_parameters) then
     print 1002
     if (nctrl == 0) then
        print 1003, eps
     else
        print 1004, nctrl
     endif
  endif
 1002 format(3x,"- Constructing B-Spline representation of surface contours")
 1003 format(8x,"multifit with automatic refinement for chisq / dof < ",e8.3)
 1004 format(8x,"multifit with ",i0," segments"/)


  ! construct B-Spline representation of surface contours
  allocate (this%slice(0:nphi-1))
  min_nctrl = huge(1)
  max_nctrl = 0
  chisq_dof(-1) = huge(1.d0)
  chisq_dof(1) = 0.d0
  if (print_progress_bar) call progress_bar(0, nphi)
  do i=0,nphi-1
     if (this%section(i)%points%nelements() < npoints) then
        call this%section(i)%savetxt("ERROR_FLUXSURF3D_POINTS")
        call ERROR("flux surface does not appear to be closed")
     endif

     if (nctrl == 0) then
        this%slice(i) = this%section(i)%bspline_autofit(k, eps, theta0, knot_balancing, fit_method)
     else
        this%slice(i) = this%section(i)%bspline_multifit(nctrl, k, theta0, knot_balancing, fit_method, eps, lambda1, lambda2)
     endif
     chisq_dof(0) = this%slice(i)%e(0) / this%slice(i)%nctrl
     min_nctrl = min(min_nctrl, this%slice(i)%nctrl)
     max_nctrl = max(max_nctrl, this%slice(i)%nctrl)
     chisq_dof(-1) = min(chisq_dof(-1), chisq_dof(0))
     chisq_dof(1) = max(chisq_dof(1), chisq_dof(0))
     if (print_progress_bar) call progress_bar(i+1, nphi)
  enddo
  if (print_parameters) then
     if (min_nctrl == max_nctrl) then
        print 1005, min_nctrl
     else
        print 1006, min_nctrl, max_nctrl
     endif
     print 1007, chisq_dof(-1), chisq_dof(1)
  endif
 1005 format(3x,"-> number of B-spline control points: ",i0)
 1006 format(3x,"-> number of B-spline control points: ",i0," - ",i0)
 1007 format(8x,"achieved chisq / dof: ",e12.4," -> ",e12.4)


  ! map parametrization
  param_ = "magnetic_angle";   if (present(param)) param_ = param
  select case (param_)
  case ("magnetic_angle")
     ! nothing to be done here

  case ("arclength")
     do i=0,nphi-1
        call this%slice(i)%set_arclength_parametrization()
     enddo

  case default
     call ERROR("invalid parametrization "//trim(param_))
  end select

  end function construct_fluxsurf3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function loadnc_fluxsurf3d(filename) result(this)
  use moose_netcdf
  use moose_utils, only: str
  use moose_geometry
  character(len=*), intent(in) :: filename
  type(fluxsurf3d)             :: this

  type(netcdf_dataset) :: nc, grp
  integer :: i, nphi


  nc = netcdf_open(filename)
  call nc%get_att("nphi", nphi)
  allocate (this%section(0:nphi-1), this%slice(0:nphi-1))
  do i=0,nphi-1
     this%section(i) = readnc_poincare_map(nc%group("poincare_map"//str(i)))
     this%slice(i)%bspline_curve = readnc_bspline_curve(nc%group("bspline"//str(i)))
  enddo
  call nc%close()

  end function loadnc_fluxsurf3d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine save(this, basename)
  class(fluxsurf3d), intent(in) :: this
  character(len=*),  intent(in) :: basename

  character(len=len_trim(basename)+8) :: filename
  integer :: i, iu


  do i=0,this%nphi-1
     ! write Poincare maps
     write (filename, 1001) trim(basename), i
     open  (newunit=iu, file=filename)
     write (iu, '(dt)') this%section(i)
     close (iu)

     ! write B-Spline interpolation
     write (filename, 1002) trim(basename), i
     call this%slice(i)%plot(filename)


     !DEBUG
     write (filename, 9001) trim(basename), i
     call this%slice(i)%P%savetxt(filename)
 9001 format(a,"_",i0,".ctr")
  enddo
 1001 format(a,"_",i0,".dat")
 1002 format(a,"_",i0,".plt")

  end subroutine save
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine savenc(this, filename)
  use moose_netcdf
  use moose_utils, only: str
  class(fluxsurf3d), intent(in) :: this
  character(len=*),  intent(in) :: filename

  type(netcdf_dataset) :: nc, grp
  integer :: i


  nc = netcdf_create(filename)
  call nc%put_att("nphi", this%nphi)
  call nc%enddef()

  do i=0,this%nphi-1
     call nc%redef()
     call nc%def_grp("poincare_map"//str(i), grp)
     call this%section(i)%writenc(grp)

     call nc%redef()
     call nc%def_grp("bspline"//str(i), grp)
     call this%slice(i)%writenc(grp)
  enddo
  call nc%close()

  end subroutine savenc
  !-----------------------------------------------------------------------------

end module flare_fluxsurf3d
