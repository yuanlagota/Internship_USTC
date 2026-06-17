module flare_m3dc1
  use iso_fortran_env
  use moose_analysis, only: poloidal_flux => scalar_mfunc2d, &
                            scalar_mfunc3d, &
                            ufunc
  use flare_bfield
  use fusion_io
  use flare_equi2d
  implicit none
  private


  public :: &
     FIO_DENSITY, &
     FIO_TEMPERATURE, &
     FIO_ELECTRON, &
     FIO_MAIN_ION


  ! M3D-C1 source file .........................................................
  type, public :: m3dc1_source
     integer :: isrc
     type(fio_search_hint), pointer :: hint

     character(len=256), private :: filename
     integer, private :: timeslice

     contains
     procedure :: broadcast => broadcast_m3dc1_source
     procedure :: free      => free_m3dc1_source
  end type m3dc1_source


  interface m3dc1_source
     procedure :: load_m3dc1_source
  end interface m3dc1_source
  ! m3dc1_src ..................................................................



  ! M3D-C1 scalar field ........................................................
  type, extends(scalar_mfunc3d), public :: m3dc1_scalar_field
     type(m3dc1_source), pointer :: src
     integer :: ifield

     logical, private :: standalone

     contains
     procedure :: free    => free_m3dc1_scalar_field

     procedure :: eval    => eval_m3dc1_scalar_field
     procedure :: deriv   => deriv_m3dc1_scalar_field
     procedure :: hessian => hessian_m3dc1_scalar_field
  end type m3dc1_scalar_field


  interface m3dc1_scalar_field
     procedure :: init_m3dc1_scalar_field
     procedure :: load_m3dc1_scalar_field
  end interface m3dc1_scalar_field
  ! m3dc1_scalar_field .........................................................



  ! M3D-C1 magnetic field ......................................................
  type, extends(magnetic_field), public :: m3dc1_bfield
     real(real64) :: phase
     integer      :: ifield

     type(m3dc1_source), private :: src
     real(real64), private :: factor

     contains
     procedure :: broadcast => broadcast_m3dc1_bfield
     procedure :: free      => free_m3dc1_bfield

     procedure :: eval
     procedure :: jac
  end type m3dc1_bfield


  interface m3dc1_bfield
     procedure :: load_m3dc1_bfield
  end interface
  ! m3dc1_bfield ...............................................................



  ! wrapper for poloidal magnetic flux .........................................
  type, extends(poloidal_flux), public :: m3dc1_psi
     type(m3dc1_source), pointer :: src
     integer :: iB, iPsi
     type(fio_search_hint), pointer :: hint

     logical, private :: standalone

     contains
     procedure :: broadcast => broadcast_m3dc1_psi
     procedure :: free => free_m3dc1_psi

     procedure :: eval => eval_m3dc1_psi
     procedure :: deriv
     procedure :: hessian
  end type m3dc1_psi


  interface m3dc1_psi
     procedure :: init_m3dc1_psi
     procedure :: load_m3dc1_psi
  end interface m3dc1_psi
  ! m3dc1_psi ..................................................................



  ! M3D-C1 equilibrium field ...................................................
  type, extends(equi2d), public :: m3dc1_equi2d
     integer :: iB, iPsi

     type(m3dc1_source), private :: src
     real(real64), private :: factor

     contains
     procedure :: broadcast => broadcast_m3dc1_equi2d
     procedure :: free      => free_m3dc1_equi2d

     procedure :: eval      => eval_m3dc1_equi2d
     procedure :: jac       => jac_m3dc1_equi2d
     procedure :: FpsiN
     procedure :: Fx
     procedure :: FdF
  end type m3dc1_equi2d
  ! m3dc1_equi2d ...............................................................


  public :: &
     load_m3dc1_equi2d

  contains
  !-----------------------------------------------------------------------------


! type m3dc1_source ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_m3dc1_source(this)
  class(m3dc1_source), intent(inout) :: this

  integer :: ierr


  ! open
  call fio_open_source_f(FIO_M3DC1_SOURCE, trim(this%filename), this%isrc, ierr)
  if (ierr /= 0) then
     print 9001, trim(this%filename)
     stop
  endif
 9001 format("ERROR: failed to open M3D-C1 source at ",a,"!")


  ! set options
  call fio_get_options_f(this%isrc, ierr)
  if (ierr /= 0) then
     print 9002;   stop
  endif
  call set_int_option(FIO_TIMESLICE, this%timeslice)
 9002 format("ERROR: failed to get M3D-C1 options!")


  ! initialize search hint
  allocate (this%hint)
  call fio_allocate_search_hint_f(this%isrc, this%hint, ierr)
  if (ierr /= 0) then
     print 9003;   stop
  endif
 9003 format("ERROR: failed to initialize search hint for M3D-C1 source!")

  end subroutine aux_init_m3dc1_source
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_m3dc1_source(filename, timeslice) result(this)
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: timeslice
  type(m3dc1_source)           :: this


  this%filename  = filename
  this%timeslice = timeslice
  call aux_init_m3dc1_source(this)

  end function load_m3dc1_source
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_m3dc1_source(this)
  use moose_mpi
  class(m3dc1_source), intent(inout) :: this

  character(len=256) :: filename
  integer :: timeslice


  call proc(0)%broadcast(this%filename)
  call proc(0)%broadcast(this%timeslice)
  if (rank > 0) call aux_init_m3dc1_source(this)

  end subroutine broadcast_m3dc1_source
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_m3dc1_source(this)
  class(m3dc1_source), intent(inout) :: this

  integer :: ierr


  call fio_deallocate_search_hint_f(this%isrc, this%hint, ierr)
  if (ierr /= 0) then
     print 9001;   stop
  endif
  deallocate (this%hint)
 9001 format("ERROR: failed to deallocate M3D-C1 search hint!")


  call fio_close_source_f(this%isrc, ierr)
  if (ierr /= 0) then
     print 9002;   stop
  endif
 9002 format("ERROR: failed to close M3D-C1 source!")

  end subroutine free_m3dc1_source
  !-----------------------------------------------------------------------------
! type m3dc1_source ============================================================



! type m3dc1_scalar_field ======================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_m3dc1_scalar_field(this, itype, ispecies)
  use moose_analysis, only: init_scalar_mfunc3d
  class(m3dc1_scalar_field), intent(inout) :: this
  integer,                   intent(in)    :: itype
  integer,                   intent(in), optional :: ispecies

  real(real64) :: r0, r1, phi0, phi1, z0, z1
  integer :: ierr


  ! set species
  if (present(ispecies)) call set_int_option(FIO_SPECIES, ispecies)


  ! get field
  this%ifield = get_field(this%src, itype)
!  call fio_field_extent_f(this%ifield, r0, r1, phi0, phi1, z0, z1, ierr)
!  if (ierr /= 0) then
!     print 9000
!     stop
!  endif
 9000 format("ERROR in m3dc1_scalar_field constructor: field extent is undefined!")
  r0 =  0.d0
  r1 =  huge(1.d0)
  z0 = -huge(1.d0)
  z1 =  huge(1.d0)
  call init_scalar_mfunc3d(this, [r0, z0, phi0], [r1, z1, phi1], [.false., .false., .true.])

  end subroutine aux_init_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function init_m3dc1_scalar_field(src, itype, ispecies) result(this)
  type(m3dc1_source), pointer, intent(in) :: src
  integer,                     intent(in) :: itype
  integer,                     intent(in), optional :: ispecies
  type(m3dc1_scalar_field)                :: this


  this%src        => src
  this%standalone = .false.
  call aux_init_m3dc1_scalar_field(this, itype, ispecies)

  end function init_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_m3dc1_scalar_field(filename, timeslice, itype, ispecies) result(this)
  character(len=*),  intent(in) :: filename
  integer,           intent(in) :: timeslice, itype
  integer,           intent(in), optional :: ispecies
  type(m3dc1_scalar_field)      :: this



  allocate (this%src, source=m3dc1_source(filename, timeslice))
  this%standalone = .true.
  call aux_init_m3dc1_scalar_field(this, itype, ispecies)

  end function load_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free_m3dc1_scalar_field(this)
  class(m3dc1_scalar_field), intent(inout) :: this


  if (this%standalone) call this%src%free()
  call close_field(this%ifield)
  call this%mfunc_free()

  end subroutine free_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_m3dc1_scalar_field(this, x) result(F)
  class(m3dc1_scalar_field), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: F

  real(real64) :: x132(3), F1(1)
  integer :: ierr


  F = 0.d0
  call MAP(x, x132)
  call fio_eval_field_f(this%ifield, x132, F1, ierr, hint=this%src%hint)
  if (ierr /= 0) return

  F = F1(1)

  end function eval_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv_m3dc1_scalar_field(this, x) result(dF)
  class(m3dc1_scalar_field), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: dF(this%ndim)

  real(real64) :: x132(3), dF132(3)
  integer :: ierr


  dF = 0.d0
  call MAP(x, x132)
  call fio_eval_field_f(this%ifield, x132, dF132, ierr, hint=this%src%hint)
  if (ierr /= 0) return

  call MAP(dF132, dF)

  end function deriv_m3dc1_scalar_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function hessian_m3dc1_scalar_field(this, x) result(H)
  class(m3dc1_scalar_field), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: H(this%ndim, this%ndim)


  ! @ todo: compute approximation of Hessian from 1st order derivatives
  H = 0.d0

  end function hessian_m3dc1_scalar_field
  !-----------------------------------------------------------------------------
! type m3dc1_scalar_field ======================================================



! type m3dc1_bfield ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_m3dc1_bfield(this)
  class(m3dc1_bfield), intent(inout) :: this


  ! set FIO_PART option to FIO_PERTURBED_ONLY
  call set_int_option(FIO_PART, FIO_PERTURBED_ONLY)


  ! apply scale factor
  if (this%factor /= 1.d0) then
     call set_real_option(FIO_LINEAR_SCALE, this%factor)
  endif


  ! read magnetic field
  this%ifield = get_field(this%src, FIO_MAGNETIC_FIELD)

  end subroutine aux_init_m3dc1_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_m3dc1_bfield(filename, timeslice, amplitude, phase) result(this)
  use moose_utils, only: basename
  use moose_math,  only: pi
  character(len=*), intent(in) :: filename
  integer,          intent(in) :: timeslice
  real(real64),     intent(in), optional :: amplitude, phase
  type(m3dc1_bfield)           :: this


  print *
  print 1000, basename(filename)
 1000 format("- Magnetic field from M3D-C1: ", a)


  call init_magnetic_field(this)
  this%factor = 1.d0;   if (present(amplitude)) this%factor = amplitude
  this%phase  = 0.d0;   if (present(phase))     this%phase  = phase
  this%src    = m3dc1_source(filename, timeslice)
  call aux_init_m3dc1_bfield(this)
  print 1001, timeslice, this%factor, this%phase
  this%phase = this%phase / 180 * pi
 1001 format(8x,"timeslice: ",i4,", amplitude = ",f0.3,", phase = ",f0.3," deg")

  end function load_m3dc1_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_m3dc1_bfield(this)
  use moose_mpi
  class(m3dc1_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  call proc(0)%broadcast(this%factor)
  call proc(0)%broadcast(this%phase)
  call this%src%broadcast()
  if (rank > 0) call aux_init_m3dc1_bfield(this)

  end subroutine broadcast_m3dc1_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_m3dc1_bfield(this)
  class(m3dc1_bfield), intent(inout) :: this

  integer :: ierr


  call close_field(this%ifield)
  call this%src%free()
  call this%mfunc_free()

  end subroutine free_m3dc1_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(m3dc1_bfield), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: B(this%mdim)

  real(real64) :: x132(3), B132(3)
  integer :: ierr


  call MAP(x, x132)
  x132(2) = x132(2) - this%phase
  call fio_eval_field_f(this%ifield, x132, B132, ierr, hint=this%src%hint)
  call MAP(B132, B)

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(m3dc1_bfield), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: jac(this%mdim, this%ndim)

  real(real64) :: x132(3), db(9)
  integer      :: ierr


  jac = 0.d0
  call MAP(x, x132)
  x132(2) = x132(2) - this%phase
  call fio_eval_field_deriv_f(this%ifield, x132, db, ierr, hint=this%src%hint)
  if (ierr /= 0) return

  jac(1,1) = db(FIO_DR_R)
  jac(2,1) = db(FIO_DR_Z)
  jac(3,1) = db(FIO_DR_PHI)
  jac(1,2) = db(FIO_DZ_R)
  jac(2,2) = db(FIO_DZ_Z)
  jac(3,2) = db(FIO_DZ_PHI)
  jac(1,3) = db(FIO_DPHI_R)
  jac(2,3) = db(FIO_DPHI_Z)
  jac(3,3) = db(FIO_DPHI_PHI)

  end function jac
  !-----------------------------------------------------------------------------
! type m3dc1_bfield ============================================================



! type m3dc1_psi ===============================================================
! constructors:
  !-----------------------------------------------------------------------------
  function init_m3dc1_psi(iB, iPsi, hint, rmin,  rmax, zmin, zmax) result(this)
  use moose_analysis, only: init_scalar_mfunc2d
  integer,                        intent(in) :: iB, iPsi
  type(fio_search_hint), pointer, intent(in) :: hint
  real(real64),                   intent(in), optional :: rmin, rmax, zmin, zmax
  type(m3dc1_psi)                            :: this


  call init_scalar_mfunc2d(this)
  nullify(this%src)
  this%iB   = iB
  this%iPsi = iPsi
  this%hint => hint
  this%standalone = .false.
  if (present(rmin)) this%lb(1) = rmin
  if (present(rmax)) this%ub(1) = rmax
  if (present(zmin)) this%lb(2) = zmin
  if (present(zmax)) this%ub(2) = zmax

  end function init_m3dc1_psi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_m3dc1_psi(filename, lb, ub) result(this)
  use moose_analysis, only: init_scalar_mfunc2d
  character(len=*), intent(in) :: filename
  real(real64),     intent(in), optional :: lb(2), ub(2)
  type(m3dc1_psi) :: this


  call init_scalar_mfunc2d(this, lb, ub)
  allocate (this%src, source=m3dc1_source(filename, 0))


  call set_int_option(FIO_PART, FIO_EQUILIBRIUM_ONLY)
  this%iPsi = get_field(this%src, FIO_VECTOR_POTENTIAL)
  this%iB   = get_field(this%src, FIO_MAGNETIC_FIELD)
  this%hint => this%src%hint
  this%standalone = .true.

  end function load_m3dc1_psi
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_m3dc1_psi(this)
  !
  ! dummy subroutine, broadcast must be done manually
  !
  class(m3dc1_psi), intent(inout) :: this
  end subroutine broadcast_m3dc1_psi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_m3dc1_psi(this)
  class(m3dc1_psi), intent(inout) :: this


  if (this%standalone) then
     call close_field(this%iB)
     call close_field(this%iPsi)
     call this%src%free()
  endif
  call this%mfunc_free()

  end subroutine free_m3dc1_psi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_m3dc1_psi(this, x) result(Psi)
  class(m3dc1_psi), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: Psi

  real(real64) :: x132(3), A(3)
  integer :: ierr


  Psi = 0.d0
  call MAP([x, 0.d0], x132)
  call fio_eval_field_f(this%iPsi, x132, A, ierr, hint=this%hint)
  if (ierr /= 0) return

  Psi = A(2) * x(1)

  end function eval_m3dc1_psi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function deriv(this, x) result(dPsi)
  class(m3dc1_psi), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: dPsi(this%ndim)

  real(real64) :: x132(3), B(3)
  integer :: ierr


  dPsi = 0.d0
  call MAP([x, 0.d0], x132)
  call fio_eval_field_f(this%iB, x132, B, ierr, hint=this%hint)
  if (ierr /= 0) return

  dPsi(1) =  B(3) * x(1)
  dPsi(2) = -B(1) * x(1)

  end function deriv
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function hessian(this, x) result(H)
  class(m3dc1_psi), intent(in) :: this
  real(real64),     intent(in) :: x(this%ndim)
  real(real64)                 :: H(this%ndim, this%ndim)

  real(real64) :: x132(3), B(3), db(9)
  integer      :: ierr


  H = 0.d0
  call MAP([x, 0.d0], x132)
  call fio_eval_field_f(this%iB, x132, B, ierr, hint=this%hint);
  if (ierr /= 0) return

  call fio_eval_field_deriv_f(this%iB, x132, db, ierr, hint=this%hint)
  if (ierr /= 0) return


  H(1,1) = B(3) + x(1) * db(FIO_DR_Z)
  H(2,2) = - x(1) * db(FIO_DZ_R)
  H(2,1) = x(1) * db(FIO_DZ_Z)
  H(1,2) = H(2,1)

  end function hessian
  !-----------------------------------------------------------------------------
! type m3dc1_psi ===============================================================



! type m3dc1_equi2d ============================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_m3dc1_equi2d(this)
  type(m3dc1_equi2d), intent(inout) :: this


  ! set FIO_PART option to FIO_EQUILIBRIUM_ONLY
  call set_int_option(FIO_PART, FIO_EQUILIBRIUM_ONLY)


  ! apply scale factor (optional)
  if (this%factor /= 1.d0) then
     call set_real_option(FIO_LINEAR_SCALE, this%factor)
  endif


  ! read vector potential and magnetic field
  this%iPsi = get_field(this%src, FIO_VECTOR_POTENTIAL)
  this%iB   = get_field(this%src, FIO_MAGNETIC_FIELD)

  end subroutine aux_init_m3dc1_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_bounding_box(this)
  use moose_math, only: pi2
  type(m3dc1_equi2d), intent(inout) :: this

  integer, parameter :: n = 360

  real(real64) :: rmin, rmax, zmin, zmax, x0(2), x1(2), theta
  integer :: i


  x0(1) = eval_series(get_series(this%src, FIO_MAGAXIS_R), 0.d0)
  x0(2) = eval_series(get_series(this%src, FIO_MAGAXIS_Z), 0.d0)
  rmin = huge(1.d0)
  rmax = 0.d0
  zmax = -rmin
  zmin = rmin
  do i=0,n-1
     theta = pi2 * i / n
     x1 = line_search(x0, x0(1) * [cos(theta), sin(theta)], 1.d-8)
     rmax = max(rmax, x1(1))
     rmin = min(rmin, x1(1))
     zmax = max(zmax, x1(2))
     zmin = min(zmin, x1(2))
  enddo
  call init_magnetic_field(this, rmin, rmax, zmin, zmax)

  contains
  !.............................................................................
  function line_search(x0, dx, eps) result(x)
  real(real64), intent(in) :: x0(2), dx(2), eps
  real(real64)             :: x(2)

  real(real64) :: x1(2), x2(2), B(3)
  integer :: i, ierr


  ! start at x0, and take coarse steps of dx to find bracket
  x2 = x0
  do
     x2 = x2 + dx
     call fio_eval_field_f(this%iB, [x2(1), 0.d0, x2(2)], B, ierr, hint=this%src%hint)
     if (ierr /= 0) exit
  enddo


  ! refinement within [x1, x2]
  x1 = x2 - dx
  do i=1,int(log(sqrt(sum(dx**2))/eps)/log(2.d0))
     x = (x2+x1) / 2
     call fio_eval_field_f(this%iB, [x(1), 0.d0, x(2)], B, ierr, hint=this%src%hint)
     if (ierr /= 0) then
        x2 = x
     else
        x1 = x
     endif
  enddo

  end function line_search
  !.............................................................................
  end subroutine aux_init_bounding_box
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_m3dc1_equi2d(filename, factor) result(this)
  !use moose_geometry, only: hypersurf2d
  use moose_utils,    only: dirname
  character(len=*),  intent(in) :: filename
  real(real64),      intent(in), optional :: factor
  type(m3dc1_equi2d)            :: this

  real(real64) :: rmin, rmax, zmin, zmax, Raxis, Zaxis, psib


  this%factor = 1.d0;   if (present(factor)) this%factor = factor
  this%src    = m3dc1_source(filename, 0)
  call aux_init_m3dc1_equi2d(this)
  call aux_init_bounding_box(this)
  rmin = this%lb(1);   rmax = this%ub(1)
  zmin = this%lb(2);   zmax = this%ub(2)
  allocate (this%Psi, source=m3dc1_psi(this%iB, this%iPsi, this%src%hint, rmin, rmax, zmin, zmax))
  print 1000, rmin, rmax, zmin, zmax
 1000 format(8x,'Bounding box:               ', &
             'R      = ',f8.3, ' m  ->  ',f8.3,' m',/36x, &
             'Z      = ',f8.3, ' m  ->  ',f8.3,' m')

  ! initialize magnetic axis
  Raxis  = eval_series(get_series(this%src, FIO_MAGAXIS_R), 0.d0)
  Zaxis  = eval_series(get_series(this%src, FIO_MAGAXIS_Z), 0.d0)
  psib   = eval_series(get_series(this%src, FIO_LCFS_PSI),  0.d0)
  this%delta_Psi = psib - this%Psi%eval([Raxis, Zaxis])
  call this%setup(Raxis, Zaxis, dirname(filename))
  print 2001, this%Bt_axis
  print 2002, this%delta_Psi
 2001 format(8x,'Toroidal magnetic field:  Bt(axis) = ',f8.3,' T')
 2002 format(8x,'Delta Psi                          = ',e12.5)

  end function load_m3dc1_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast_m3dc1_equi2d(this)
  use moose_mpi
  class(m3dc1_equi2d), intent(inout) :: this


  if (rank > 0) allocate (m3dc1_psi :: this%Psi)
  call this%equi2d_broadcast()
  call this%src%broadcast()
  call proc(0)%broadcast(this%factor)
  if (rank > 0) then
     call aux_init_m3dc1_equi2d(this)
     this%Psi = m3dc1_psi(this%iB, this%iPsi, this%src%hint)
     this%Psi%lb = this%lb(1:2)
     this%Psi%ub = this%ub(1:2)
  endif

  end subroutine broadcast_m3dc1_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free_m3dc1_equi2d(this)
  class(m3dc1_equi2d), intent(inout) :: this


  call close_field(this%iB)
  call close_field(this%iPsi)
  call this%src%free()
  call this%mfunc_free()

  end subroutine free_m3dc1_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_m3dc1_equi2d(this, x) result(B)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),              intent(in) :: x(this%ndim)
  real(real64)                          :: B(this%mdim)

  real(real64) :: x132(3), B132(3)
  integer :: ierr


  B = 0.d0
  call MAP(x, x132)
  call fio_eval_field_f(this%iB, x132, B132, ierr, hint=this%src%hint)
  if (ierr /= 0) return

  call MAP(B132, B)

  end function eval_m3dc1_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac_m3dc1_equi2d(this, x) result(jac)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: x(this%ndim)
  real(real64)                    :: jac(this%mdim, this%ndim)

  real(real64) :: x132(3), db(9)
  integer      :: ierr


  jac = 0.d0
  call MAP(x, x132)
  call fio_eval_field_deriv_f(this%iB, x132, db, ierr, hint=this%src%hint)
  if (ierr /= 0) return

  jac(1,1) = db(FIO_DR_R)
  jac(2,1) = db(FIO_DR_Z)
  jac(3,1) = db(FIO_DR_PHI)
  jac(1,2) = db(FIO_DZ_R)
  jac(2,2) = db(FIO_DZ_Z)
  jac(3,2) = db(FIO_DZ_PHI)
  jac(1,3) = db(FIO_DPHI_R)
  jac(2,3) = db(FIO_DPHI_Z)
  jac(3,3) = db(FIO_DPHI_PHI)

  end function jac_m3dc1_equi2d
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function Fx(this, r) result(F)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: F

  real(real64) :: B(3)


  B = this%eval(r)
  F = B(3) * r(1)

  end function Fx
  !-----------------------------------------------------------------------------
  function FpsiN(this, psiN) result(F)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: psiN
  real(real64)                    :: F

  real(real64) :: B(3)


  print *, "ERROR: FpsiN not implemented for M3D-C1 equilibrium"
  stop

  end function FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function FdF(this, r)
  class(m3dc1_equi2d), intent(in) :: this
  real(real64),        intent(in) :: r(:)
  real(real64)                    :: FdF(0:1)

  real(real64) :: B(3), J(3,3)


  B      = this%eval(r)
  J      = this%jac(r)
  FdF(0) = B(3) * r(1)
  FdF(1) = - J(3,2) * r(1)**2 * B(1)

  end function FdF
  !-----------------------------------------------------------------------------
! type m3dc1_equi2d ============================================================



! module procedures:
  !-----------------------------------------------------------------------------
  subroutine close_field(ifield)
  integer, intent(in) :: ifield

  integer :: ierr


  call fio_close_field_f(ifield, ierr)
  if (ierr /= 0) then
     print 9000;   stop
  endif
 9000 format("ERROR: failed to close M3D-C1 field!")

  end subroutine close_field
  !-----------------------------------------------------------------------------
  function get_field(src, ifield_type) result(ifield)
  type(m3dc1_source), intent(in) :: src
  integer,            intent(in) :: ifield_type
  integer                        :: ifield

  integer :: ierr


  call fio_get_field_f(src%isrc, ifield_type, ifield, ierr)
  if (ierr /= 0) then
     print 9000, ifield_type;   stop
  endif
 9000 format("ERROR: failed to get M3D-C1 field ",i0,"!")

  end function get_field
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine set_int_option(iopt, val)
  integer,          intent(in) :: iopt, val

  integer :: ierr


  call fio_set_int_option_f(iopt, val, ierr)
  if (ierr /= 0) then
     print 9000, iopt;   stop
  endif
 9000 format("ERROR: failed to set M3D-C1 option ",i0,"!")

  end subroutine set_int_option
  !-----------------------------------------------------------------------------
  subroutine set_real_option(iopt, val)
  integer,          intent(in) :: iopt
  real(real64),     intent(in) :: val

  integer :: ierr


  call fio_set_real_option_f(iopt, val, ierr)
  if (ierr /= 0) then
     print 9000, iopt;   stop
  endif
 9000 format("ERROR: failed to set M3D-C1 option ",i0,"!")

  end subroutine set_real_option
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval_series(iseries, t) result(s)
  integer,      intent(in) :: iseries
  real(real64), intent(in) :: t
  real(real64)             :: s

  integer :: ierr


  call fio_eval_series_f(iseries, t, s, ierr)
  if (ierr /= 0) then
     print 9000, iseries;   stop
  endif
 9000 format("ERROR: failed to eval M3D-C1 series ",i0,"!")

  end function eval_series
  !-----------------------------------------------------------------------------
  function get_series(src, iseries_type) result(iseries)
  type(m3dc1_source), intent(in) :: src
  integer,            intent(in) :: iseries_type
  integer                        :: iseries

  integer :: ierr


  call fio_get_series_f(src%isrc, iseries_type, iseries, ierr)
  if (ierr /= 0) then
     print 9000, iseries_type;   stop
  endif
 9000 format("ERROR: failed to get M3D-C1 series ",i0,"!")

  end function get_series
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine MAP(xin, xout)
  real(real64), intent(in)  :: xin(3)
  real(real64), intent(out) :: xout(3)


  xout(1) = xin(1)
  xout(2) = xin(3)
  xout(3) = xin(2)

  end subroutine MAP
  !-----------------------------------------------------------------------------

end module flare_m3dc1
