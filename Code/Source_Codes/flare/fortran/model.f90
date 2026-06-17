!===============================================================================
! Definition of magnetic field and boundary geometry.
!===============================================================================
module flare_model
  use iso_fortran_env
  use moose_geometry, only: hypersurf3d, hypersurf2d
  use flare_control,  only: report
  use flare_bfield
  use flare_bspline3d
  use flare_coilset
  use flare_harmonic
  use flare_hint
  use flare_interp
  use flare_jorek
  use flare_m3dc1
  use flare_marsf
  use flare_equi2d,   only: equi2d_bfield => equi2d
  implicit none
  private


  character(len=*), parameter, public :: &
     TYPE_EQUI2D_AMHD      = "amhd", &
     TYPE_EQUI2D_GEQDSK    = "geqdsk", &
     TYPE_EQUI2D_SONNET    = "sonnet", &
     TYPE_EQUI2D_JOREK     = "jorek", &
     TYPE_EQUI2D_M3DC1     = "m3dc1", &
     TYPE_EQUI3D           = "equi3d", &
     TYPE_EQUI3D_BMW       = "equi3d_bmw", &
     TYPE_EQUI3D_MGRID     = "equi3d_mgrid", &
     TYPE_BFIELD_BSPLINE3D = "bspline3d", &
     TYPE_BFIELD_COILSET   = "coilset", &
     TYPE_BFIELD_GPEC      = "gpec", &
     TYPE_BFIELD_INTERP    = "interp", &
     TYPE_BFIELD_JOREK     = "jorek", &
     TYPE_BFIELD_M3DC1     = "m3dc1", &
     TYPE_BFIELD_MARSF     = "marsf"



  type :: bfield_container
     class(magnetic_field), allocatable :: implementation
  end type bfield_container



  ! magnetic field model (equilibrium + perturbation fields) ...................
  type, extends(magnetic_field), public :: flare_magnetic_field
     ! equilibrium magnetic field
     class(equilibrium_bfield), allocatable :: equi
     character(len=128) :: type_equi2d = ""

     ! perturbation field
     type(bfield_container), allocatable :: perturbation(:)
     integer   :: nperturbation = 0

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: eval, perturbation_eval
     procedure :: jac
  end type flare_magnetic_field
  ! flare_magnetic_field .......................................................



  ! model data
  type(flare_magnetic_field), target, public :: bfield
  type(hypersurf3d), target, public :: boundary
  class(equi2d_bfield), pointer, public :: equi2d => null()
  type(hypersurf2d), pointer, public :: boundary2d => null()



  public :: &
     assert_model, assert_equi2d, &
     alloc_boundary, load_axisurf, load_torosurf, loadnc_boundary, &
     alloc_equi2d, load_geqdsk, load_equi2d_jorek, load_equi2d_m3dc1, load_sonnet, load_amhd, &
     alloc_equi3d, load_equi3d_bmw, load_equi3d_coilset, load_equi3d_hint, load_equi3d_interp, load_equi3d_mgrid, &
     alloc_perturbation, load_bspline3d, load_coilset, load_gpec, load_interp, load_jorek, load_m3dc1, load_marsf, &
     load_model, broadcast_model, setup_model, free_model

  contains
  !-----------------------------------------------------------------------------


! type flare_magnetic_field ====================================================
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  use flare_amhd
  use flare_equi2d, only: geqdsk_equi2d, sonnet_equi2d
  use flare_equi3d
  class(flare_magnetic_field), intent(inout) :: this

  character(len=32) :: bfield_type
  integer :: i, k = 0


  call this%bfield_broadcast()


  ! equilibrium field
  call proc(0)%broadcast(this%type_equi2d)
  if (rank > 0) then
     select case(this%type_equi2d)
     case(TYPE_EQUI2D_AMHD)
        allocate (amhd_equi2d   :: this%equi)

     case(TYPE_EQUI2D_GEQDSK)
        allocate (geqdsk_equi2d :: this%equi)

     case(TYPE_EQUI2D_SONNET)
        allocate (sonnet_equi2d :: this%equi)

     case(TYPE_EQUI2D_JOREK)
        allocate (jorek_equi2d  :: this%equi)

     case(TYPE_EQUI2D_M3DC1)
        allocate (m3dc1_equi2d  :: this%equi)

     case(TYPE_EQUI3D)
        allocate (equi3d        :: this%equi)

     end select
  endif
  call this%equi%broadcast()


  ! perturbation field
  call proc(0)%broadcast(this%nperturbation)
  if (rank > 0  .and.  this%nperturbation > 0) allocate (this%perturbation(this%nperturbation))
  do i=1,this%nperturbation
     associate (bfield => this%perturbation(i)%implementation)
     if (rank == 0) then
        select type(bfield)
        type is(bspline3d_bfield)
           bfield_type = TYPE_BFIELD_BSPLINE3D
        type is(coilset)
           bfield_type = TYPE_BFIELD_COILSET
        type is(gpec_bfield)
           bfield_type = TYPE_BFIELD_GPEC
        type is(interp_bfield)
           bfield_type = TYPE_BFIELD_INTERP
        type is(jorek_bfield)
           bfield_type = TYPE_BFIELD_JOREK
           k = bfield%k
        type is(m3dc1_bfield)
           bfield_type = TYPE_BFIELD_M3DC1
        type is(marsf_bfield)
           bfield_type = TYPE_BFIELD_MARSF
        end select
     endif
     end associate
     call proc(0)%broadcast(bfield_type)
     call proc(0)%broadcast(k)

     if (rank > 0) then
        select case(bfield_type)
        case(TYPE_BFIELD_BSPLINE3D)
           allocate (bspline3d_bfield :: this%perturbation(i)%implementation)
        case(TYPE_BFIELD_COILSET)
           allocate (coilset          :: this%perturbation(i)%implementation)
        case(TYPE_BFIELD_GPEC)
           allocate (gpec_bfield      :: this%perturbation(i)%implementation)
        case(TYPE_BFIELD_INTERP)
           allocate (interp_bfield    :: this%perturbation(i)%implementation)
        case(TYPE_BFIELD_JOREK)
           allocate (jorek_bfield     :: this%perturbation(i)%implementation)
           select type(P => this%perturbation(i)%implementation)
           type is(jorek_bfield)
              select type(E => this%equi)
              type is(jorek_equi2d)
                 call P%aux_broadcast(E, k)
              end select
           end select

        case(TYPE_BFIELD_M3DC1)
           allocate (m3dc1_bfield     :: this%perturbation(i)%implementation)
        case(TYPE_BFIELD_MARSF)
           allocate (marsf_bfield     :: this%perturbation(i)%implementation)
        end select
     endif
     call this%perturbation(i)%implementation%broadcast()
  enddo

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(flare_magnetic_field), intent(inout) :: this

  integer :: i


  ! equilibrium field
  call this%equi%free()
  deallocate (this%equi)


  ! perturbation field
  do i=1,this%nperturbation
     call this%perturbation(i)%implementation%free()
  enddo
  if (this%nperturbation > 0) deallocate (this%perturbation)


  ! finalize
  call this%mfunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  use flare_control, only: BFIELD_EVAL, increment_counter
  class(flare_magnetic_field), intent(in) :: this
  real(real64),                intent(in) :: x(this%ndim)
  real(real64)                            :: B(this%mdim)


  B = this%equi%eval(x) + this%perturbation_eval(x)
  call increment_counter(BFIELD_EVAL)

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function perturbation_eval(this, x) result(B)
  class(flare_magnetic_field), intent(in) :: this
  real(real64),                intent(in) :: x(this%ndim)
  real(real64)                            :: B(this%mdim)

  integer :: i


  B = 0.d0
  do i=1,this%nperturbation
     B = B + this%perturbation(i)%implementation%eval(x)
  enddo

  end function perturbation_eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  use flare_control, only: JACOBIAN_EVAL, increment_counter
  class(flare_magnetic_field), intent(in) :: this
  real(real64),                intent(in) :: x(this%ndim)
  real(real64)                            :: jac(this%mdim, this%ndim)

  integer :: i


  jac = this%equi%jac(x)
  do i=1,this%nperturbation
     jac = jac + this%perturbation(i)%implementation%jac(x)
  enddo
  call increment_counter(JACOBIAN_EVAL)

  end function jac
  !-----------------------------------------------------------------------------
! type flare_magnetic_field ====================================================



! module procedures:
  !-----------------------------------------------------------------------------
  subroutine assert_model()
  use moose_error


  if (.not.allocated(bfield%equi)) call ERROR("model is not defined")

  end subroutine assert_model
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine assert_equi2d(procedure_name)
  use moose_error
  use moose_mpi
  character(len=*), intent(in) :: procedure_name


  call assert_model()
  select type(equi => bfield%equi)
  class is(equi2d_bfield)
     equi2d => equi
     if (.not.associated(boundary2d)) allocate (boundary2d, source=boundary%rzslice(0.d0))

  class default
     if (rank == 0) call ERROR("equilibrium is not axisymmetric", procedure_name)
  end select

  end subroutine assert_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_equi2d(type_equi2d, source)
  character(len=*),     intent(in) :: type_equi2d
  class(equi2d_bfield), intent(in) :: source


  allocate (bfield%equi, source = source)
  bfield%type_equi2d = type_equi2d

  end subroutine aux_init_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_equi3d(i, source)
  use moose_error
  use flare_equi3d
  integer,               intent(in) :: i
  class(magnetic_field), intent(in) :: source


  select type(equi => bfield%equi)
  type is (equi3d)
     if (i <= 0  .or.  i > equi%n) then
        print *, "i, n = ", i, equi%n
        call ERROR("equilibrium field out of range")
     endif
     allocate (equi%container(i)%implementation, source=source)
  end select

  end subroutine aux_init_equi3d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_perturbation(i, source)
  use moose_error
  integer,               intent(in) :: i
  class(magnetic_field), intent(in) :: source


  if (i <= 0  .or.  i > bfield%nperturbation) then
     print *, "i, n = ", i, bfield%nperturbation
     call ERROR("perturbation field out of range")
  endif
  allocate (bfield%perturbation(i)%implementation, source=source)

  end subroutine aux_init_perturbation
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function get_model_path(model, database) result(path)
  use moose_error
  use moose_utils,   only: expanduser, join, startswith, isdir
  use flare_control, only: database_dict => database
  character(len=*), intent(in) :: model
  character(len=*), intent(in), optional :: database
  character(len=:), allocatable :: path


  ! current working directory or explicit path
  if (model == ""  .or.  startswith(model, "/")  .or.  startswith(model, "./")) then
     path = model

  ! path relative to database
  elseif (present(database)) then
     if (.not.database_dict%has_key(database)) then
        call ERROR("database '"//trim(database)//"' is not defined")
     endif
     path = join(expanduser(database_dict%get(database)), model)

  else
     path = join(expanduser(database_dict%get("default")), model)
  endif


  ! check if path exists
  if (path /= "") then
     if (.not.isdir(path)) call ERROR("model "//trim(model)//" does not exist")
  endif

  end function get_model_path
  !-----------------------------------------------------------------------------


! boundary =====================================================================
  !-----------------------------------------------------------------------------
  subroutine alloc_boundary(nexpl, nimpl)
  !
  ! nexpl:   number of surfaces defined through configuration file
  ! nimpl:   number of surfaces through stellarator symmetry
  !
  use moose_hypersurface, only: alloc_hypersurf3d
  integer, intent(in) :: nexpl, nimpl


  if (report) then
     print 1000
     print *, "Boundary setup (divertor targets, first wall, limiters):"
     print *
  endif
 1000 format(80("="))

  boundary = alloc_hypersurf3d(nexpl + nimpl)
  boundary%nsurfaces = nexpl   ! reset nsurfaces - will be incremented in load_torosurf

  end subroutine alloc_boundary
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_axisurf(i, key, filename, units)
  use moose_axisurf
  integer, intent(in) :: i
  character(len=*), intent(in) :: key, filename, units


  boundary%surfaces(i)%key = key
  allocate (boundary%surfaces(i)%geometry, source = axisurf(filename, units, convert_units="m"))
  select type (S => boundary%surfaces(i)%geometry)
  type is (axisurf)
     if (report) then
        print 2001, i, key, S%description(), S%P%segments() + 1
        print *
     endif
  end select
 2001 format(3x,"- ",i0,": ",a," - ",a,/,8x,"Number of nodes: ",i0)

  end subroutine load_axisurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine load_torosurf(i, key, filename, stellarator_symmetry, units, vfallback)
  use moose_math
  use moose_torosurf
  integer,          intent(in) :: i
  character(len=*), intent(in) :: key, filename, units
  logical,          intent(in) :: stellarator_symmetry
  character(len=*), intent(in), optional :: vfallback

  real(real64), parameter :: eps_phi = 1.d-9

  real(real64) :: phi1
  integer :: isym


  boundary%surfaces(i)%key = key
  allocate (boundary%surfaces(i)%geometry, source = torosurf(filename, units, convert_units="m", vfallback=vfallback))
  select type (S => boundary%surfaces(i)%geometry)
  type is (torosurf)
     ! snap to lower bound of field period
     if (abs(S%phi(0)) < eps_phi) S%phi(0) = 0.d0
     phi1 = pi2 / S%symmetry

     if (stellarator_symmetry) then
        ! snap to upper bound of half-field period
        if (abs(S%phi(S%nu) - phi1/2) < eps_phi) S%phi(S%nu) = phi1/2

        ! add stellarator symmetric surface
        boundary%nsurfaces = boundary%nsurfaces + 1
        isym = boundary%nsurfaces
        boundary%surfaces(isym)%key = boundary%surfaces(i)%key // ".T"
        allocate (boundary%surfaces(isym)%geometry, source = stellarator_symmetric_torosurf(S))

     else
        ! snap to upper bound of field period
        if (abs(S%phi(S%nu) - phi1) < eps_phi) S%phi(S%nu) = phi1
     endif


     if (report) then
        if (stellarator_symmetry) then
           print 2001, i, isym, key, S%description(), S%nu + 1, S%nv + 1
        else
           print 2002, i, key, S%description(), S%nu + 1, S%nv + 1
        endif
        print *
     endif
  end select
 2001 format(3x,"- ",i0,",",i0,": ",a," - ",a,/,8x,"Number of nodes: ",i0," x ",i0)
 2002 format(3x,"- ",i0,": ",a," - ",a,/,8x,"Number of nodes: ",i0," x ",i0)

  contains
  !.............................................................................
  function stellarator_symmetric_torosurf(S) result(this)
  type(torosurf), intent(in) :: S
  type(torosurf)             :: this

  integer :: j, k


  this = torosurf(S%nu, S%nv, S%symmetry)
  call this%metadata%set("description", S%description() // "(stellarator symmetric)")
  call this%metadata%set("units", S%units())
  do j=0,S%nu
     this%phi(j) = phi1 - S%phi(S%nu-j)
     do k=0,S%nv
        this%rz(1, k, j) =  S%rz(1, S%nv-k, S%nu-j)
        this%rz(2, k, j) = -S%rz(2, S%nv-k, S%nu-j)
        this%v(k, j) = -S%v(S%nv-k, S%nu-j)
     enddo
  enddo
  this%phi_order = -S%phi_order
  call setup_torosurf(this)


  ! snap to bounds of half-field period
  if (abs(this%phi(this%nu) - phi1) < eps_phi) this%phi(this%nu) = phi1
  if (abs(this%phi(0) - phi1/2) < eps_phi) this%phi(0) = phi1/2

  end function stellarator_symmetric_torosurf
  !.............................................................................
  end subroutine load_torosurf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine loadnc_boundary(filename)
  use moose_geometry
  character(len=*), intent(in) :: filename

  integer :: i


  boundary = loadnc_hypersurf3d(filename, convert_units="m")
  if (report) then
     print 1000
     print *, "Boundary setup (divertor targets, first wall, limiters):"
     print *

     do i=1,boundary%nsurfaces
        select type (S => boundary%surfaces(i)%geometry)
        type is (axisurf)
           print 2001, boundary%surfaces(i)%key, S%description(), S%P%segments() + 1

        type is (torosurf)
           print 2002, boundary%surfaces(i)%key, S%description(), S%nu + 1, S%nv + 1
        end select
        print *
     enddo
  endif
 1000 format(80("="))
 2001 format(3x,"- ",i0,": ",a," - ",a,/,8x,"Number of nodes: ",i0)
 2002 format(3x,"- ",i0,": ",a," - ",a,/,8x,"Number of nodes: ",i0," x ",i0)

  end subroutine loadnc_boundary
  !-----------------------------------------------------------------------------
! boundary =====================================================================



! magnetic field ===============================================================
  !-----------------------------------------------------------------------------
  subroutine alloc_equi2d()


  if (report) then
     print *
     print 1000
     print *, "Toroidally symmetric (2D) equilibrium:"
  endif
 1000 format(80("="))

  end subroutine alloc_equi2d
  !-----------------------------------------------------------------------------
  subroutine load_geqdsk(filename, scale_Bt, scale_Ip, spline_order)
  use flare_equi2d, only: load_geqdsk_equi2d
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale_Bt, scale_Ip
  integer,          intent(in), optional :: spline_order


  call aux_init_equi2d(TYPE_EQUI2D_GEQDSK, load_geqdsk_equi2d(filename, scale_Bt, scale_Ip, spline_order))

  end subroutine load_geqdsk
  !-----------------------------------------------------------------------------
  subroutine load_equi2d_jorek(filename)
  character(len=*), intent(in) :: filename


  call aux_init_equi2d(TYPE_EQUI2D_JOREK, load_jorek_equi2d(filename))

  end subroutine load_equi2d_jorek
  !-----------------------------------------------------------------------------
  subroutine load_equi2d_m3dc1(filename)
  character(len=*), intent(in) :: filename


  call aux_init_equi2d(TYPE_EQUI2D_M3DC1, load_m3dc1_equi2d(filename))

  end subroutine load_equi2d_m3dc1
  !-----------------------------------------------------------------------------
  subroutine load_sonnet(filename, scale_Bt, scale_Ip)
  use flare_equi2d, only: load_sonnet_equi2d
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: scale_Bt, scale_Ip


  call aux_init_equi2d(TYPE_EQUI2D_SONNET, load_sonnet_equi2d(filename, scale_Bt, scale_Ip))

  end subroutine load_sonnet
  !-----------------------------------------------------------------------------
  subroutine load_amhd(R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX)
  use flare_amhd
  real(real64), intent(in) :: R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX


  call aux_init_equi2d(TYPE_EQUI2D_AMHD, init_amhd_equi2d(R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX))

  end subroutine load_amhd
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine alloc_equi3d(nequi3d, axis3d)
  use moose_utils, only: dirname
  use flare_equi3d
  integer,          intent(in) :: nequi3d
  character(len=*), intent(in) :: axis3d


  if (report) then
     print *
     print 1000
     print *, "Non-axisymmetric (3D) equilibrium:"
  endif
 1000 format(80("="))

  allocate (bfield%equi, source=new_equi3d(nequi3d, axis3d, dirname(axis3d)//"/.psiN"))
  bfield%type_equi2d = TYPE_EQUI3D

  end subroutine alloc_equi3d
  !-----------------------------------------------------------------------------
  subroutine load_equi3d_bmw(i, filename, amplitude, spline_order)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  integer,          intent(in), optional :: spline_order


  call aux_init_equi3d(i, load_bmw_bfield(filename, amplitude, spline_order))

  end subroutine load_equi3d_bmw
  !-----------------------------------------------------------------------------
  subroutine load_equi3d_coilset(i, filename, amplitude, units)
  use moose_r3grid, only: length_scale
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename, units
  real(real64),     intent(in) :: amplitude


  call aux_init_equi3d(i, coilset(filename, amplitude, length_scale(units)))

  end subroutine load_equi3d_coilset
  !-----------------------------------------------------------------------------
  subroutine load_equi3d_hint(i, filename, group, bmax)
  integer,          intent(in) :: i, group
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: bmax


  call aux_init_equi3d(i, load_hint_bfield(filename, group, bmax))

  end subroutine load_equi3d_hint
  !-----------------------------------------------------------------------------
  subroutine load_equi3d_interp(i, filename, filetype, amplitude, bfield_units, length_units)
  use moose_r3grid, only: length_scale
  use flare_bfield, only: bfield_scale
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename, filetype, bfield_units, length_units
  real(real64),     intent(in) :: amplitude

  real(real64) :: bscale, lscale

  bscale = bfield_scale(bfield_units)
  lscale = length_scale(length_units)
  call aux_init_equi3d(i, load_interp_bfield(filename, filetype, amplitude, bscale, lscale))

  end subroutine load_equi3d_interp
  !-----------------------------------------------------------------------------
  subroutine load_equi3d_mgrid(i, filename, amplitudes, dtype, spline_order)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitudes(:)
  character(len=*), intent(in), optional :: dtype
  integer,          intent(in), optional :: spline_order


  call aux_init_equi3d(i, load_mgrid_bfield(filename, amplitudes, dtype, spline_order))

  end subroutine load_equi3d_mgrid
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine alloc_perturbation(n)
  integer, intent(in) :: n


  bfield%nperturbation = n
  if (n <= 0) return
  allocate (bfield%perturbation(n))
  if (report) then
     print *
     print 1000
     print *, "Perturbation field:"
  endif
 1000 format(80("="))

  end subroutine alloc_perturbation
  !-----------------------------------------------------------------------------
  subroutine load_bspline3d(i, filename, amplitude, dtype, spline_order, value_order)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  character(len=*), intent(in), optional :: dtype, value_order
  integer,          intent(in), optional :: spline_order


  call aux_init_perturbation(i, load_bspline3d_bfield(filename, amplitude, dtype, spline_order, value_order))

  end subroutine load_bspline3d
  !-----------------------------------------------------------------------------
  subroutine load_coilset(i, filename, amplitude, units)
  use moose_r3grid, only: length_scale
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename, units
  real(real64),     intent(in) :: amplitude


  call aux_init_perturbation(i, coilset(filename, amplitude, length_scale(units)))

  end subroutine load_coilset
  !-----------------------------------------------------------------------------
  subroutine load_gpec(i, filename, amplitude, phase)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  real(real64),     intent(in), optional :: phase


  call aux_init_perturbation(i, load_gpec_bfield(filename, amplitude, phase))

  end subroutine load_gpec
  !-----------------------------------------------------------------------------
  subroutine load_interp(i, filename, filetype, amplitude, bfield_units, length_units)
  use moose_r3grid, only: length_scale
  use flare_bfield, only: bfield_scale
  integer,          intent(in) :: i
  character(len=*), intent(in) :: filename, filetype, bfield_units, length_units
  real(real64),     intent(in) :: amplitude

  real(real64) :: bscale, lscale


  bscale = bfield_scale(bfield_units)
  lscale = length_scale(length_units)
  call aux_init_perturbation(i, load_interp_bfield(filename, filetype, amplitude, bscale, lscale))

  end subroutine load_interp
  !-----------------------------------------------------------------------------
  subroutine load_jorek(i, n, amplitude)
  use moose_error
  integer,          intent(in) :: i, n
  real(real64),     intent(in) :: amplitude


  select type(E => bfield%equi)
  class is(jorek_equi2d)
     call aux_init_perturbation(i, init_jorek_bfield(n, amplitude, E))

  class default
     call ERROR("Perturbation field from JOREK is incompatible with equilibrium")
  end select

  end subroutine load_jorek
  !-----------------------------------------------------------------------------
  subroutine load_m3dc1(i, filename, timeslice, amplitude, phase)
  integer,          intent(in) :: i, timeslice
  character(len=*), intent(in) :: filename
  real(real64),     intent(in) :: amplitude
  real(real64),     intent(in), optional :: phase


  call aux_init_perturbation(i, m3dc1_bfield(filename, timeslice, amplitude, phase))

  end subroutine load_m3dc1
  !-----------------------------------------------------------------------------
  subroutine load_marsf(i, schimesh, bplasma, amplitude, phase)
  integer,          intent(in) :: i
  character(len=*), intent(in) :: schimesh, bplasma
  real(real64),     intent(in) :: amplitude
  real(real64),     intent(in), optional :: phase


  call aux_init_perturbation(i, marsf_bfield(schimesh, bplasma, amplitude, phase))

  end subroutine load_marsf
  !-----------------------------------------------------------------------------
! magnetic field ===============================================================


  !-----------------------------------------------------------------------------
  subroutine load_model(model, database)
  !
  ! load selected *model* from *database*
  !
  use moose_error
  use moose_configparser
  use moose_utils, only: join, isdir, startswith
  character(len=*), intent(in) :: model
  character(len=*), intent(in), optional :: database

  integer, parameter :: &
     ALLOC_STAGE = 1, &
     LOAD_STAGE = 2

  type(configparser) :: config
  type(cp_section), pointer :: section
  character(len=:), allocatable :: model_path, config_path, key, dtype
  integer :: ipass, nboundary, nZ, nequi2d, nequi3d, nperturbation


  model_path = get_model_path(model, database)


  ! 1. boundary geometry
  call read_config(model_path, ".boundary")
  do ipass=ALLOC_STAGE,LOAD_STAGE
     nboundary = 0
     nZ = 0
     section => config%first_section()
     do
        if (.not.associated(section)) exit
        call get_dtype(section%key, key, dtype)
        call load_boundary(key, dtype)

        section => config%next_section()
     enddo

     ! allocate memory
     if (ipass == ALLOC_STAGE) call alloc_boundary(nboundary, nZ)
  enddo


  ! 2. magnetic field definition
  call read_config(model_path, ".bfield")
  ! 2.1. equilibrium
  do ipass=ALLOC_STAGE,LOAD_STAGE
     nequi2d = 0
     nequi3d = 0
     section => config%first_section()
     do
        if (.not.associated(section)) exit
        call get_dtype(section%key, key, dtype)
        ! axisymmetric (2D) equilibrium
        if (startswith(dtype, "equi2d_")) then
           nequi2d = nequi2d + 1
           if (ipass == LOAD_STAGE) call load_equi2d(dtype(8:))

        ! non-axisymmetric (3D) equilibrium
        elseif (startswith(dtype, "equi3d_")) then
           nequi3d = nequi3d + 1
           if (ipass == LOAD_STAGE) call load_equi3d(dtype)
        endif

        section => config%next_section()
     enddo

     ! allocate memory
     if (ipass == ALLOC_STAGE) then
        if (nequi2d > 1) call ERROR("multiple definitions of toroidally symmetric equilibria not allowed")
        if (nequi2d == 1) then
           call alloc_equi2d()
        else
           call alloc_equi3d(nequi3d, join(config_path, ".equi3d"))
        endif
     endif

  enddo

  ! 2.2. perturbation
  do ipass=ALLOC_STAGE,LOAD_STAGE
     nperturbation = 0
     section => config%first_section()
     do
        if (.not.associated(section)) exit
        call get_dtype(section%key, key, dtype)
        nperturbation = nperturbation + 1
        if (ipass == LOAD_STAGE) call load_perturbation(dtype)

        section => config%next_section()
     enddo

     ! allocate memory
     if (ipass == ALLOC_STAGE) call alloc_perturbation(nperturbation)
  enddo
  call setup_model()

  contains
  !.............................................................................
  subroutine read_config(model_path, suffix)
  use moose_error
  use moose_utils, only: isdir, join
  character(len=*), intent(in) :: model_path, suffix

  character(len=:), allocatable :: filename
  logical :: ex


  config_path = model_path
  filename = join(model_path, suffix)
  if (isdir(filename)) then
     config_path = filename
     filename = join(config_path, suffix)
  endif

  inquire(file=filename, exist=ex)
  if (.not.ex) call ERROR("configuration file "//trim(filename)//" does not exist")
  config = configparser()
  call config%read(filename)

  end subroutine read_config
  !.............................................................................
  subroutine get_dtype(section, key, dtype)
  use moose_utils, only: split
  character(len=*),              intent(in   ) :: section
  character(len=:), allocatable, intent(  out) :: key, dtype

  integer :: i, j


  i = scan(section, ':')
  if (i == 0) then
     key = trim(section)
     dtype = key
     return
  endif
  key = section(1:i-1)


  j = scan(section(i+1:), ':')
  if (j == 0) then
     dtype = trim(section(i+1:))
  else
     dtype = section(i+1:j)
  endif

  end subroutine get_dtype
  !.............................................................................
  subroutine load_boundary(key, dtype)
  character(len=*), intent(in) :: key, dtype

  character(len=:), allocatable :: filename, units
  logical :: stellarator_symmetry


  nboundary = nboundary + 1
  select case(dtype)
  case("axisurf")
     if (ipass == LOAD_STAGE) then
        filename = join(config_path, section%get("filename"))
        units = section%get("units", fallback="m")
        call load_axisurf(nboundary, key, filename, units)
     endif


  case("torosurf")
     stellarator_symmetry = section%getlogical("stellarator_symmetry", fallback=.false.)
     if (stellarator_symmetry) nZ = nZ + 1
     if (ipass == LOAD_STAGE) then
        filename = join(config_path, section%get("filename"))
        units = section%get("units", fallback="m")
        call load_torosurf(nboundary, key, filename, stellarator_symmetry, units)
     endif


  case default
     call ERROR("invalid boundary type '"//dtype//"'")
  end select

  end subroutine load_boundary
  !.............................................................................
  subroutine load_equi2d(equi2d_dtype)
  character(len=*), intent(in) :: equi2d_dtype

  character(len=:), allocatable :: filename
  real(real64) :: scale_Bt, scale_Ip, R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX
  integer :: spline_order


  select case(equi2d_dtype)
  case (TYPE_EQUI2D_AMHD)
     R0           = section%getdouble("R0")
     Z0           = section%getdouble("Z0")
     Bt           = section%getdouble("Bt")
     Ip           = section%getdouble("Ip")
     A            = section%getdouble("A")
     eps          = section%getdouble("eps")
     kappa        = section%getdouble("kappa")
     delta        = section%getdouble("delta")
     rX           = section%getdouble("rX")
     zX           = section%getdouble("zX")
     call load_amhd(R0, Z0, Bt, Ip, A, eps, kappa, delta, rX, zX)


  case (TYPE_EQUI2D_GEQDSK)
     filename     = join(config_path, section%get("filename"))
     scale_Bt     = section%getdouble("scale_Bt",     fallback = 1.d0)
     scale_Ip     = section%getdouble("scale_Ip",     fallback = 1.d0)
     spline_order = section%getint   ("spline_order", fallback = 4)
     call load_geqdsk(filename, scale_Bt, scale_Ip, spline_order)


  case (TYPE_EQUI2D_JOREK)
     filename     = join(config_path, section%get("filename"))
     call load_equi2d_jorek(filename)


  case (TYPE_EQUI2D_M3DC1)
     filename     = join(config_path, section%get("filename"))
     call load_equi2d_m3dc1(filename)


  case (TYPE_EQUI2D_SONNET)
     filename     = join(config_path, section%get("filename"))
     scale_Bt     = section%getdouble("scale_Bt",     fallback = 1.d0)
     scale_Ip     = section%getdouble("scale_Ip",     fallback = 1.d0)
     call load_sonnet(filename, scale_Bt, scale_Ip)


  case default
     call ERROR("invalid equi2d type '"//equi2d_dtype//"'")
  end select
  call section%remove();  deallocate (section)

  end subroutine load_equi2d
  !.............................................................................
  subroutine load_equi3d(equi3d_dtype)
  use flare_equi3d
  character(len=*), intent(in) :: equi3d_dtype

  character(len=:), allocatable :: bfield_units, dtype, filename, filetype, length_units, units
  real(real64), allocatable :: amplitudes(:)
  real(real64) :: amplitude, bmax
  integer :: group, spline_order


  select case(equi3d_dtype)
  case (TYPE_EQUI3D_BMW)
     filename     = join(config_path, section%get("filename"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     spline_order = section%getint   ("spline_order", fallback = 5)
     call load_equi3d_bmw(nequi3d, filename, amplitude, spline_order)


  case (TYPE_EQUI3D_COILSET)
     filename     = join(config_path, section%get("filename"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     units        = section%get      ("units",        fallback = "m")
     call load_equi3d_coilset(nequi3d, filename, amplitude, units)


  case (TYPE_EQUI3D_HINT)
     filename     = join(config_path, section%get("filename"))
     group        = section%getint   ("group",        fallback = -1)
     bmax         = section%getdouble("bmax",         fallback = 0.d0)
     call load_equi3d_hint(nequi3d, filename, group, bmax)


  case (TYPE_EQUI3D_INTERP)
     filename     = join(config_path, section%get("filename"))
     filetype     = section%get      ("filetype",     fallback = "ascii")
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     bfield_units = section%get      ("bfield_units", fallback = "T")
     length_units = section%get      ("length_units", fallback = "m")
     call load_equi3d_interp(nequi3d, filename, filetype, amplitude, bfield_units, length_units)


  case (TYPE_EQUI3D_MGRID)
     filename     = join(config_path, section%get("filename"))
     amplitudes   = section%getarray ("amplitudes")
     dtype        = section%get      ("dtype",        fallback = "magnetic_field")
     spline_order = section%getint   ("spline_order", fallback = 0)
     call load_equi3d_mgrid(nequi3d, filename, amplitudes, dtype, spline_order)


  case default
     call ERROR("invalid equi3d type '"//equi3d_dtype//"'")
  end select
  call section%remove();  deallocate (section)

  end subroutine load_equi3d
  !.............................................................................
  subroutine load_perturbation(perturbation_dtype)
  character(len=*), intent(in) :: perturbation_dtype

  character(len=:), allocatable :: bfield_units, bplasma, dtype, filename, &
     filetype, length_units, schimesh, units, value_order
  real(real64) :: amplitude, phase
  integer :: n, spline_order, timeslice


  select case(perturbation_dtype)
  case (TYPE_BFIELD_BSPLINE3D)
     filename     = join(config_path, section%get("filename"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     dtype        = section%get      ("dtype",        fallback = "vector_potential")
     spline_order = section%getint   ("spline_order", fallback = 0)
     value_order  = section%get      ("value_order",  fallback = "column_major")
     call load_bspline3d(nperturbation, filename, amplitude, dtype, spline_order, value_order)


  case (TYPE_BFIELD_COILSET)
     filename     = join(config_path, section%get("filename"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     units        = section%get      ("units",        fallback = "m")
     call load_coilset(nperturbation, filename, amplitude, units)


  case (TYPE_BFIELD_GPEC)
     filename     = join(config_path, section%get("filename"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     phase        = section%getdouble("phase",        fallback = 0.d0)
     call load_gpec(nperturbation, filename, amplitude, phase)


  case (TYPE_BFIELD_INTERP)
     filename     = join(config_path, section%get("filename"))
     filetype     = section%get      ("filetype",     fallback = "ascii")
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     bfield_units = section%get      ("bfield_units", fallback = "T")
     length_units = section%get      ("length_units", fallback = "m")
     call load_interp(nperturbation, filename, filetype, amplitude, bfield_units, length_units)


  case (TYPE_BFIELD_JOREK)
     n            = section%getint   ("n")
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     call load_jorek(nperturbation, n, amplitude)


  case (TYPE_BFIELD_M3DC1)
     filename     = join(config_path, section%get("filename"))
     timeslice    = section%getint   ("timeslice",    fallback = 0)
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     phase        = section%getdouble("phase",        fallback = 0.d0)
     call load_m3dc1(nperturbation, filename, timeslice, amplitude, phase)


  case (TYPE_BFIELD_MARSF)
     schimesh     = join(config_path, section%get("schimesh"))
     bplasma      = join(config_path, section%get("bplasma"))
     amplitude    = section%getdouble("amplitude",    fallback = 1.d0)
     phase        = section%getdouble("phase",        fallback = 0.d0)
     call load_marsf(nperturbation, schimesh, bplasma, amplitude, phase)

  case default
     call ERROR("invalid bfield type '"//perturbation_dtype//"'")
  end select

  end subroutine load_perturbation
  !.............................................................................
  end subroutine load_model
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine broadcast_model()


  call boundary%broadcast()
  call bfield%broadcast()

  end subroutine broadcast_model
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine setup_model()
  use moose_error
  use moose_math, only: gcd
  use flare_equi3d
  real(real64) :: lb(2), ub(2)
  integer :: i, k, nfp


  if (bfield%type_equi2d == "") call ERROR("equilibrium is undefined")
  associate (equi => bfield%equi)
  select type(equi)
  type is (equi3d)
     call equi%setup()
  end select
  end associate


  ! set bounding box and symmetry
  lb = bfield%equi%lb(1:2)
  ub = bfield%equi%ub(1:2)
  nfp = bfield%equi%nfp
  do i=1,bfield%nperturbation
     do k=1,2
        lb(k) = max(lb(k), bfield%perturbation(i)%implementation%lb(k))
        ub(k) = min(ub(k), bfield%perturbation(i)%implementation%ub(k))
     enddo
     nfp = gcd(nfp, bfield%perturbation(i)%implementation%nfp)
  enddo
  call init_magnetic_field(bfield, lb(1), ub(1), lb(2), ub(2), nfp)
  if (report) then
     print *
     print 1000
  endif
 1000 format(80("="))

  end subroutine setup_model
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_model()


  call boundary%free()
  call bfield%free()
  if (associated(boundary2d)) then
     call boundary2d%free()
     deallocate (boundary2d)
  endif

  end subroutine free_model
  !-----------------------------------------------------------------------------

end module flare_model
