!===============================================================================
! Axisymmetric equilibrium fields (tokamaks)
!===============================================================================
module flare_equi2d
  use iso_fortran_env
  use moose_rlist
  use moose_analysis, only: scalar_mfunc2d, interp
  use flare_control,  only: report
  use flare_bfield,   only: equilibrium_bfield
  implicit none
  private


! plain data:
  ! geqdsk dataset .............................................................
  type, public :: geqdsk
     character(len=48) :: label
     integer :: nr, nz, nbbbs, limitr
     real(real64) :: Rdim, Zdim, Rcentr, Rleft, Zmid, Rmaxis, Zmaxis, Simag, &
                     Sibry, Bcentr, Current

     real(real64), allocatable :: fpol(:), pres(:), ffprim(:), pprime(:), psirz(:,:), &
                                  qpsi(:), rbbbs(:), zbbbs(:), rlim(:), zlim(:)
  end type geqdsk
  ! geqdsk .....................................................................



! bfield implementations:
  ! base class for axisymmetric equilibrium fields .............................
  type, abstract, extends(equilibrium_bfield), public :: equi2d
     ! poloidal flux function [Wb/rad]
     class(scalar_mfunc2d), allocatable :: Psi

     ! position of magnetic axis [m]
     real(real64) :: r0(2)
     ! toroidal magnetic field and poloidal flux on axis, and difference to boundary
     real(real64) :: Bt_axis, Psi_axis, delta_Psi

     ! X-points
     type(rlist2) :: X

     contains
     procedure :: find_xpoints
     procedure :: equi2d_broadcast
     procedure :: equi2d_free
     procedure :: free => equi2d_free
     procedure :: setup, aux_setup_magnetic_axis

     procedure :: magnetic_axis
     procedure :: psiN, grad_psiN
     procedure :: psiN_X
     generic   :: F => Fx, FpsiN
     procedure :: Fx
     procedure(FpsiN), deferred :: FpsiN
     procedure(FdF), deferred :: FdF
     procedure :: eval
     procedure :: jac
     procedure :: dpsi

     procedure :: xpoint, xpoint_hessian, xpoint_stability
     procedure :: kappa
     procedure :: rzcoords, rzcoordsX
  end type equi2d


  abstract interface
     function FpsiN(this, psiN) result(F)
     !
     ! toroidal field function [T m]
     !
     import equi2d, real64
     class(equi2d), intent(in) :: this
     real(real64),  intent(in) :: psiN
     real(real64)              :: F
     end function FpsiN


     function FdF(this, r)
     !
     ! toroidal field function [T m] and its first derivative with respect to psiN
     !
     import equi2d, real64
     class(equi2d), intent(in) :: this
     real(real64),  intent(in) :: r(:)
     real(real64)              :: FdF(0:1)
     end function FdF
  end interface
  ! equi2d .....................................................................



  ! geqdsk implementation of equi2d ............................................
  type, extends(equi2d), public :: geqdsk_equi2d
     ! poloidal current function [T m]
     type(interp) :: Ffunc

     ! internal values for poloidal flux on boundary and axis (-> eval Ffunc)
     real(real64) :: Sibry, Simag

     contains
     procedure :: broadcast => geqdsk_equi2d_broadcast
     procedure :: free      => geqdsk_equi2d_free

     procedure :: Fx        => geqdsk_equi2d_Fx
     procedure :: FpsiN     => geqdsk_equi2d_FpsiN
     procedure :: FdF       => geqdsk_equi2d_FdF
  end type geqdsk_equi2d
  ! geqdsk_equi2d ..............................................................



  ! sonnet implementation of equi2d ............................................
  type, extends(equi2d), public :: sonnet_equi2d
     real(real64) :: F0

     contains
     procedure :: broadcast => sonnet_equi2d_broadcast

     procedure :: FpsiN     => sonnet_equi2d_FpsiN
     procedure :: FdF       => sonnet_equi2d_FdF
  end type sonnet_equi2d
  ! sonnet_equi2d ..............................................................



  public :: &
     load_geqdsk, &
     load_geqdsk_equi2d, &
     load_sonnet_equi2d

  contains
  !-----------------------------------------------------------------------------


! type geqdsk ==================================================================
  !-----------------------------------------------------------------------------
  function load_geqdsk(filename) result(this)
  !
  ! load plain geqdsk data file
  !
  use ieee_arithmetic
  use moose_error
  character(len=*), intent(in) :: filename
  type(geqdsk)                 :: this

  character(len=80) :: s
  character(len=8)  :: r5fmt
  real(real64)      :: dummy
  integer           :: i, iu, ierr, j, nr, nz


  open  (newunit=iu, file=filename, action="read")
  ! read header
  read  (iu, '(a)') s
  this%label = s(:48);   read  (s(49:), *, iostat=ierr) i, this%nr, this%nz
  if (ierr /= 0) call ERROR("unexpected header format in geqdsk file")
  r5fmt = '(5e16.9)'
  read  (iu, '(a)') s
  do i=1,2
     read (s, r5fmt, iostat=ierr) this%Rdim, this%Zdim, this%Rcentr, this%Rleft, this%Zmid
     if (ierr == 0) exit
     r5fmt = '(5e17.9)'
  enddo
  read  (iu, r5fmt) this%Rmaxis,  this%Zmaxis,  this%Simag,   this%Sibry,   this%Bcentr
  read  (iu, r5fmt) this%Current, this%Simag,   dummy,        this%Rmaxis,  dummy
  read  (iu, r5fmt) this%Zmaxis,  dummy,        this%Sibry,   dummy,        dummy

  ! read equilibrium data
  nr = this%nr
  nz = this%nz
  allocate (this%fpol(nr), this%pres(nr), this%ffprim(nr), this%pprime(nr), &
            this%psirz(nr,nz), this%qpsi(nr))
  read  (iu, r5fmt) (this%fpol(i), i=1,nr)
  read  (iu, r5fmt) (this%pres(i), i=1,nr)
  read  (iu, r5fmt) (this%ffprim(i), i=1,nr)
  read  (iu, r5fmt) (this%pprime(i), i=1,nr)
  read  (iu, r5fmt) ((this%psirz(i,j), i=1,nr), j=1,nz)
  read  (iu, r5fmt) (this%qpsi(i), i=1,nr)

  ! device and plasma boundary
  read  (iu, *) this%nbbbs, this%limitr
  allocate (this%rbbbs(this%nbbbs), this%zbbbs(this%nbbbs), &
            this%rlim(this%limitr), this%zlim(this%limitr))
  read  (iu, r5fmt) (this%rbbbs(i), this%zbbbs(i), i=1,this%nbbbs)
  read  (iu, r5fmt) (this%rlim(i),  this%zlim(i),  i=1,this%limitr)
  close (iu)

  ! sanity check
  if (any(ieee_is_nan(this%fpol))) call ERROR("invalid fpol", "load_geqdsk")
  if (any(ieee_is_nan(this%pres))) call ERROR("invalid pres", "load_geqdsk")
  if (any(ieee_is_nan(this%ffprim))) call ERROR("invalid ffprim", "load_geqdsk")
  if (any(ieee_is_nan(this%pprime))) call ERROR("invalid pprime", "load_geqdsk")
  if (any(ieee_is_nan(this%psirz))) call ERROR("invalid psirz", "load_geqdsk")
  if (any(ieee_is_nan(this%qpsi))) call ERROR("invalid qpsi", "load_geqdsk")
  if (any(ieee_is_nan(this%rbbbs))) call ERROR("invalid rbbbs", "load_geqdsk")
  if (any(ieee_is_nan(this%zbbbs))) call ERROR("invalid zbbbs", "load_geqdsk")
  if (any(ieee_is_nan(this%rlim))) call ERROR("invalid rlim", "load_geqdsk")
  if (any(ieee_is_nan(this%zlim))) call ERROR("invalid zlim", "load_geqdsk")

  end function load_geqdsk
  !-----------------------------------------------------------------------------
! type geqdsk ==================================================================



! class equi2d =================================================================
! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_setup_magnetic_axis(this, R0, Z0)
  !
  ! initialize magnetic axis from approximate position at (R0, Z0)
  !
  ! requires:
  !    this%Psi   (poloidal magnetic flux function)
  !
  ! output:
  !    this%r0
  !    this%Psi_axis
  !
  use moose_mfunc
  class(equi2d), intent(inout) :: this
  real(real64),  intent(in)    :: R0, Z0

  integer :: istat


  call find_nearest_critical_point(this%Psi, (/R0, Z0/), this%r0, istat)
  if (istat /= 0) then
     print *, "istat = ", istat
     print 9000;   stop
  endif
 9000 format("ERROR: aux_setup_magnetic_axis failed!")


  this%Psi_axis = this%Psi%eval(this%r0)

  end subroutine aux_setup_magnetic_axis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_xpoints(this, dirname)
  !
  ! load X-point setup (if available)
  !
  class(equi2d),     intent(inout) :: this
  character(len=*),  intent(in)    :: dirname

  character(len=len_trim(dirname)+8) :: filename
  real(real64) :: x(2)
  logical :: ex
  integer :: i, iu, nx


  this%X   = rlist2()
  filename = trim(dirname)//"/.equi2d"
  inquire (file=filename, exist=ex)
  if (ex) then
     open  (newunit=iu, file=filename)
     read  (iu, *) nx
     if (nx > 0  .and.  report) print 3100
     do i=1,nx
        read  (iu, *) x
        call this%X%append(x)
        if (report) print 3101, x
     enddo
     if (report) print *
     close (iu)
  else
     if (report) print 3200
  endif
 3100 format(8x,"X-point(s):")
 3101 format(16x,"(",f6.3,", ",f6.3,") m")
 3200 format(8x,"no X-points defined")

  end subroutine aux_init_xpoints
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_init_helicity(this, stdout)
  !
  ! input:
  !    this%r0           (position of magnetic axis)
  !    sign(delta_Psi)   (direction of poloidal magnetic flux)
  !
  ! output:
  !    this%Bt_axis
  !    this%Bt_sign
  !    this%Bp_sign
  !    this%helicity
  !
  class(equi2d), intent(inout) :: this
  logical,       intent(in)    :: stdout

  character(len=1), parameter :: helicity(-1:1) = (/'-', '?', '+'/)

  real(real64) :: B(3), r0(3)


  r0(1:2) = this%r0
  r0(3)   = 0.d0
  B       = this%eval(r0)

  this%Bt_axis   = B(3)
  this%Bt_sign   = -1;   if (B(3) > 0)           this%Bt_sign =  1
  this%Bp_sign   =  1;   if (this%delta_Psi > 0) this%Bp_sign = -1
  this%helicity  = this%Bt_sign * this%Bp_sign
  if (stdout) print 1000, helicity(this%Bt_sign), helicity(this%Bp_sign)
 1000 format(8x,"Helicity (toroidal & poloidal field direction):",4x,a,2x,a)

  end subroutine aux_init_helicity
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine setup(this, Raxis0, Zaxis0, dirname)
  !
  class(equi2d),     intent(inout) :: this
  real(real64),      intent(in)    :: Raxis0, Zaxis0
  character(len=*),  intent(in)    :: dirname


  call this%aux_setup_magnetic_axis(Raxis0, Zaxis0)
  if (report) print 1000, this%r0
 1000 format(8x,"Magnetic axis:",19x,"(",f6.3,",",2x,f6.3,") m")


  call aux_init_xpoints(this, dirname)
  if (this%X%nelements() > 0) this%delta_Psi = this%Psi%eval(this%X%element(0)) - this%Psi_axis
  call aux_init_helicity(this, report)

  end subroutine setup
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine find_xpoints(this, boundary, dirname, nsample, rrange, zrange)
  !
  ! initialize X-points
  !
  ! input:
  !    this%Psi
  !
  ! output:
  !    this%X
  !
  use moose_utils, only: str
  use moose_linalg
  use moose_hypersurface
  use moose_mfunc, only: find_critical_points2d
  use moose_contours
  use flare_control, only: step_size => separatrix2d_step_size, epsabs => separatrix2d_epsabs, &
                           offset => separatrix2d_offset, fX => separatrix2d_fX, &
                           alpha => separatrix2d_alpha, nmax => separatrix2d_nmax
  class(equi2d),     intent(inout) :: this
  type(hypersurf2d), intent(in)    :: boundary
  character(len=*),  intent(in)    :: dirname
  integer,           intent(in), optional :: nsample(2)
  real(real64),      intent(in), optional :: rrange(2), zrange(2)

  character(len=*), parameter :: case_symbol(4) = ['O', 'O', 'X', '?']

  character(len=len_trim(dirname)+8) :: filename
  type(xcontour), allocatable :: S(:)
  type(rlist2) :: Xlist
  real(real64) :: x(2), H(2,2), xrange(2,2)
  integer      :: i, icase, iu


  ! 0. user defined scan domain
  xrange(1,:) = this%lb(1:2)
  xrange(2,:) = this%ub(1:2)
  if (present(rrange)) then
     if (rrange(2) > rrange(1)) xrange(:,1) = rrange
  endif
  if (present(zrange)) then
     if (zrange(2) > zrange(1)) xrange(:,2) = zrange
  endif


  call this%X%clear()
  ! 1. scan domain for X-points
  print 1000
  Xlist  = find_critical_points2d(this%Psi, nsample, xrange)
  ! sort out minima and maxima
  i = 0
  do
     if (i == Xlist%nelements()) exit

     x = Xlist%element(i)
     H = this%Psi%hessian(x)

     icase = second_partial_derivative_test_2d(H)
     print 1001, x, case_symbol(icase)
     if (case_symbol(icase) == "X") then
        i = i + 1
     else
        call Xlist%drop(i)
     endif
  enddo
  print 1002
  if (Xlist%nelements() == 0) then
     print *, "no critical points founds"
     return
  endif
 1000 format(8x,"Scanning for critical points ...")
 1001 format(16x,"(",f6.3,", ",f6.3,") m",4x,a)
 1002 format(8x,"...done")


  ! 2. sort out irrelevant X-points (outside boundary)
  i = 0
  do
     if (i == Xlist%nelements()) exit

     if (boundary%includes(Xlist%element(i))) then
        i = i + 1
     else
        call Xlist%drop(i)
     endif
  enddo


  ! 3.a. generate separatrix for each X-point and find main one
  allocate (S(Xlist%nelements()))
  print 3000
  do i=1,Xlist%nelements()
     print 3001, i, Xlist%element(i-1)
     S(i) = xcontour(this%Psi, Xlist, i-1, step_size, boundary=boundary, epsabs=epsabs, &
               offset=offset, fX=fX, nmax=nmax, alpha=alpha)
     call S(i)%savetxt(trim(dirname)//"/separatrix"//str(i))

     if (any(S(i)%iconnect > 0)) then
        call this%X%append(S(i)%x)
        this%delta_Psi = this%Psi%eval(S(i)%x) - this%Psi_axis
     endif
  enddo
  print 1002
 3000 format(8x,"Constructing separatrices for relevant X-points ...")
 3001 format(16x,i0,": (",f6.3,", ",f6.3,") m")

  ! 3.b. add secondary X-points
  do i=1,Xlist%nelements()
     if (.not.any(S(i)%iconnect > 0)) then
        call this%X%append(S(i)%x)
     endif
  enddo

  ! 3.c. save X-point setup
  print 3100
  filename = trim(dirname)//"/.equi2d"
  open  (newunit=iu, file=filename)
  write (iu, *) this%X%nelements()
  do i=0,this%X%nelements()-1
     x = this%X%element(i)
     write (iu, *) x
     print 3101, x
  enddo
  print *
  close (iu)
 3100 format(8x,"X-point(s):")
 3101 format(16x,"(",f6.3,", ",f6.3,") m")

  end subroutine find_xpoints
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine equi2d_broadcast(this)
  use moose_mpi
  class(equi2d), intent(inout) :: this


  call this%equilibrium_broadcast()
  call this%Psi%broadcast()
  call proc(0)%broadcast(this%r0)
  call proc(0)%broadcast(this%Bt_axis)
  call proc(0)%broadcast(this%Psi_axis)
  call proc(0)%broadcast(this%delta_Psi)
  call this%X%broadcast()

  end subroutine equi2d_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine equi2d_free(this)
  use moose_mpi
  class(equi2d), intent(inout) :: this


  call this%X%free()
  call this%Psi%free()
  deallocate (this%Psi)
  call this%mfunc_free()

  end subroutine equi2d_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function magnetic_axis(this, phi) result(r0)
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: phi
  real(real64)              :: r0(2)


  r0 = this%r0

  end function magnetic_axis
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function psiN(this, r)
  !
  ! return normalized poloidal flux at r = (R[m], Z[m])
  !
  use moose_error
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: r(:)
  real(real64)              :: psiN


  if (size(r) < 2) call ERROR("size(r) < 2 not allowed", "equi2d%psiN")
  psiN = (this%Psi%eval(r) - this%Psi_axis) / this%delta_Psi

  end function psiN
  !-----------------------------------------------------------------------------
  function psiN_X(this, ix) result(psiN)
  !
  ! normalized poloidal flux at ix-th X-point
  !
  class(equi2d), intent(in) :: this
  integer,       intent(in) :: ix
  real(real64)              :: psiN


  psiN = this%psiN(this%X%element(ix-1))

  end function psiN_X
  !-----------------------------------------------------------------------------
  function grad_psiN(this, r)
  !
  ! gradient of PsiN
  !
  use moose_error
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: r(:)
  real(real64)              :: grad_psiN(3)


  if (size(r) < 2) call ERROR("size(r) < 2 not allowed", "equi2d%psiN")
  grad_psiN(1:2) = (this%Psi%deriv(r) - this%Psi_axis) / this%delta_Psi
  grad_psiN(3) = 0.d0

  end function grad_psiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function Fx(this, r)
  !
  ! return toroidal field function [T m] at r = (R[m], Z[m])
  !
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: r(:)
  real(real64)              :: Fx

  real(real64) :: psiN


  psiN = (this%Psi%eval(r(1:2)) - this%Psi_axis) / this%delta_Psi
  Fx   = this%F(psiN)

  end function Fx
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(Bf)
  !
  ! return (Br, Bz, Bphi) [T]
  !
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: x(this%ndim)
  real(real64)              :: Bf(this%ndim)

  real(real64) :: dPsidr(2)


  ! poloidal field
  dPsidr = this%Psi%deriv(x(1:2))
  Bf(1)  = -dPsidr(2) / x(1)
  Bf(2)  =  dPsidr(1) / x(1)

  ! toroidal field
  Bf(3)  = this%F(x) / x(1)

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: x(this%ndim)
  real(real64)              :: jac(this%ndim, this%ndim)

  real(real64) :: psiN, dPsidr(2), H(2,2), FdF(0:1)


  jac       = 0.d0
  psiN      = this%psiN(x(1:2));   if (psiN > 1.d0) psiN = 1.d0
  dPsidr    = this%Psi%deriv(x(1:2))
  H         = this%Psi%Hessian(x(1:2))
  FdF       = this%FdF(x)

  ! Br
  jac(1,1)  = -H(2,1)/x(1) + dPsidr(2)/x(1)**2
  jac(1,2)  = -H(2,2)/x(1)
  ! Bz
  jac(2,1)  =  H(1,1)/x(1) - dPsidr(1)/x(1)**2
  jac(2,2)  =  H(1,2)/x(1)
  ! Bphi
  jac(3,1)  =  FdF(1) * dPsidr(1) / x(1) - FdF(0)/x(1)**2
  jac(3,2)  =  FdF(1) * dPsidr(2) / x(1)

  end function jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function dpsi(this, x)
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: x(2)
  real(real64)              :: dpsi

  real(real64) :: gradPsi(2)


  gradPsi = this%Psi%deriv(x(1:2))
  dpsi    = sqrt(sum(gradPsi**2))

  end function dpsi
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xpoint(this, i)
  !
  ! return coordinates of i-th X-point
  !
  use moose_error
  use moose_utils, only: ordinal
  class(equi2d), intent(in) :: this
  integer,       intent(in) :: i
  real(real64)              :: xpoint(2)


  if (i <= 0) call ERROR("X-point index <= 0 not allowed")
  if (i > this%X%nelements()) call ERROR(ordinal(i)//" X-point is not defined, check equilibrium configuration")

  xpoint = this%X%element(i-1)

  end function xpoint
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xpoint_hessian(this, i) result(v)
  !
  ! compute eigenvalues and eigenvectors of Hessian(psiN) at i-th X-point.
  !
  use moose_linalg, only: hessian2d_analysis
  class(equi2d), intent(in) :: this
  integer,       intent(in) :: i
  real(real64)              :: v(2,-1:1)

  real(real64) :: H(2,2), lambda1, lambda2, theta, x0(2)


  x0 = this%xpoint(i)
  H  = this%Psi%hessian(x0) / this%delta_psi
  theta = this%poloidal_angle(x0)
  call hessian2d_analysis(H(1,1), H(1,2), H(2,2), lambda1, lambda2, v(:,-1), v(:,1), theta)
  v(:,0) = [lambda1, lambda2]

  end function xpoint_hessian
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function xpoint_stability(this, i) result(v)
  !
  ! compute eigenvalues and eigenvectors of poloidal field Jacobian at i-th X-point.
  !
  ! **Result:**
  !
  !   v(:,-1)    Direction of unstable manifold
  !   v(:, 0)    Eigenvalues [lambda1, lambda2]
  !   v(:, 1)    Direction of stable manifold
  !
  use moose_linalg, only: stability_analysis
  class(equi2d), intent(in) :: this
  integer,       intent(in) :: i
  real(real64)              :: v(2,-1:1)

  real(real64) :: H(2,2), lambda1, lambda2, theta, x0(2)


  x0 = this%xpoint(i)
  H  = this%Psi%hessian(x0)
  theta = this%poloidal_angle(x0)
  call stability_analysis(H, lambda1, lambda2, v(:,-1), v(:,1), theta)
  v(:,0) = [lambda1, lambda2]

  end function xpoint_stability
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function kappa(this, r)
  !
  ! return local curvature of flux surface contour
  !
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: r(2)
  real(real64)              :: kappa

  real(real64) :: dPsi(2), H(2,2), dd(2)


  dPsi  = this%Psi%deriv(r)
  H     = this%Psi%hessian(r)
  dd(1) = dPsi(2)**2 * H(1,1)  -  2 * dPsi(1) * dPsi(2) * H(1,2)  +  dPsi(1)**2 * H(2,2)
  dd(2) = (dPsi(1)**2  +  dPsi(2)**2)**(1.5d0)
  kappa = - dd(1) / dd(2)

  end function kappa
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rzcoords(this, psiN, theta, x0, ierr) result(r)
  !
  ! convert magnetic coordinates (psiN, theta[rad]) to cylindrical coordinates (r[m], z[m])
  !
  ! optional parameters:
  !    x0    initial guess for (r,z)
  !
  use moose_error
  use moose_math, only: pi
  use moose_mfunc
  class(equi2d), intent(in   ) :: this
  real(real64),  intent(in   ) :: psiN, theta
  real(real64),  intent(in   ), optional :: x0(2)
  integer,       intent(  out), optional :: ierr
  real(real64)                 :: r(2)

  character(len=256) :: msg
  real(real64) :: v(2), s, psi0
  integer :: istat


  v(1) = cos(theta)
  v(2) = sin(theta)
  ! user defined initial point
  if (present(x0)) then
     s = sum((x0 - this%r0)*v)

  ! set initial point near axis
  else
     s = 0.1d0 * this%r0(1)
  endif
  r = this%r0 + s * v


  ! find r with Psi(r) = Psi0 along direction v
  psi0 = psiN * this%delta_Psi + this%Psi_axis
  call root_finder(this%Psi, r, psi0, v, istat)
  ! re-try with stiffer iteration
  if (istat /= 0) then
     r = this%r0 + 0.1d0 * this%r0(1) * v
     call root_finder(this%Psi, r, psi0, v, istat, damping=0.5d0)
  endif

  if (present(ierr)) then
     ierr = istat
  elseif (istat /= 0) then
     write (msg, 9000) psiN, theta/pi*180.d0
     call ERROR(msg)
  endif
 9000 format("ERROR: evaluation of R-Z coordinates failed for (psiN, theta [deg]) = (",f0.3,", ",f0.3,")")

  end function rzcoords
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function rzcoordsX(this, psiN) result(rz)
  !
  ! compute (R,Z) coordinates for psiN in [0,1] on line between magnetic axis and X-point.
  !
  use moose_error
  use moose_utils, only: str
  use moose_mfunc, only: line_search
  class(equi2d), intent(in) :: this
  real(real64),  intent(in) :: psiN
  real(real64)              :: rz(2)

  real(real64) :: psiT
  integer      :: istat


  psiT = psiN * this%delta_Psi + this%Psi_axis
  rz = line_search(this%Psi, this%r0, this%xpoint(1), psiT, istat)
  if (istat /= 0) then
     call ERROR("line search failed with istat = "//str(istat))
  endif

  end function rzcoordsX
  !-----------------------------------------------------------------------------
! class equi2d =================================================================



! type geqgsk_equi2d ===========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function load_geqdsk_equi2d(filename, scale_Bt, scale_Ip, spline_order) result(this)
  use moose_math,     only: linspace
  use moose_analysis, only: bspline2d, cspline
  use moose_utils,    only: dirname
  use flare_bfield,   only: init_magnetic_field
  character(len=*), intent(in) :: filename
  type(geqdsk_equi2d)          :: this
  real(real64),     intent(in), optional :: scale_Bt, scale_Ip
  integer,          intent(in), optional :: spline_order

  real(real64), allocatable :: R(:), Z(:), Psi(:)
  type(geqdsk) :: g
  integer :: k


  g = load_geqdsk(filename)
  if (present(scale_Bt)) then
     g%Bcentr = g%Bcentr * scale_Bt
     g%fpol   = g%fpol   * scale_Bt
  endif
  if (present(scale_Ip)) then
     g%Simag  = g%Simag  * scale_Ip
     g%Sibry  = g%Sibry  * scale_Ip
     g%psirz  = g%psirz  * scale_Ip
  endif


  ! set up spatial discretization
  allocate (R(g%nr),   source=linspace(g%Rleft,            g%Rleft+g%Rdim,     g%nr))
  allocate (Z(g%nz),   source=linspace(g%Zmid-g%Zdim/2.d0, g%Zmid+g%Zdim/2.d0, g%nz))
  allocate (Psi(g%nr), source=linspace(0.d0,               1.d0,               g%nr))

  if (report) then
     print 1000, adjustl(g%label)
     print *
     print 1001, g%nr, g%nz
     print 1002, R(1), R(g%nr), Z(1), Z(g%nz)
     print 1003, g%Rcentr
     print 1004, g%Bcentr
  endif
 1000 format(8x,a)
 1001 format(8x,'Grid resolution: ',10x,i4,' x ',i4, ' nodes')
 1002 format(8x,'Computational Box:          ', &
             'R      = ',f8.3, ' m  ->  ',f8.3,' m',/36x, &
             'Z      = ',f8.3, ' m  ->  ',f8.3,' m')
 1003 format(8x,'Reference position:         R0     = ',f8.3,' m')
 1004 format(8x,'Toroidal magnetic field:    Bt(R0) = ',f8.3,' T')


  ! initialize equi2d
  call init_magnetic_field(this, R(1), R(g%nr), Z(1), Z(g%nz))
  this%Ffunc = cspline(Psi, g%fpol)
  this%Simag = g%Simag
  this%Sibry = g%Sibry
  k = 4;   if (present(spline_order)) k = spline_order
  allocate (this%Psi, source=bspline2d(R, Z, g%psirz, spline_order=k+1))

  this%delta_Psi = g%Sibry - g%Simag
  call setup(this, g%Rmaxis, g%Zmaxis, dirname(filename))
  if (report  .and.  this%Bp_sign * g%Current < 0.d0) print 8001
 8001 format(1x,"WARNING: sign of current metadata is inconsistent with poloidal flux")

  end function load_geqdsk_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine geqdsk_equi2d_broadcast(this)
  use moose_mpi
  use moose_analysis, only: bspline2d
  class(geqdsk_equi2d), intent(inout) :: this


  if (rank > 0) allocate (bspline2d :: this%Psi)
  call this%equi2d_broadcast()
  call this%Ffunc%broadcast()
  call proc(0)%broadcast(this%Simag)
  call proc(0)%broadcast(this%Sibry)

  end subroutine geqdsk_equi2d_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine geqdsk_equi2d_free(this)
  use moose_mpi
  class(geqdsk_equi2d), intent(inout) :: this


  call this%Ffunc%free()
  call this%equi2d_free()

  end subroutine geqdsk_equi2d_free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function geqdsk_equi2d_Fx(this, r) result(F)
  !
  ! return toroidal field function [T m] at r = (R[m], Z[m])
  ! NOTE: this overrides equi2d%Fx for explict use of Simag and Sibry provided by the equilibrium
  !       instead of Psi_axis and delta_Psi evaluated from the interpolated data.
  !
  class(geqdsk_equi2d), intent(in) :: this
  real(real64),         intent(in) :: r(:)
  real(real64)                     :: F

  real(real64) :: psiN


  psiN = (this%Psi%eval(r(1:2)) - this%Simag) / (this%Sibry - this%Simag)
  if (psiN > 1.d0  .or.  psiN < 0.d0) psiN = 1.d0
  F    = this%Ffunc%eval(psiN)

  end function geqdsk_equi2d_Fx
  !-----------------------------------------------------------------------------
  function geqdsk_equi2d_FpsiN(this, psiN) result(F)
  !
  ! return toroidal field function [T m]
  !
  class(geqdsk_equi2d), intent(in) :: this
  real(real64),         intent(in) :: psiN
  real(real64)                     :: F

  real(real64) :: psiN_


  psiN_ = psiN;   if (psiN > 1.d0  .or.  psiN < 0.d0) psiN_ = 1.d0
  F = this%Ffunc%eval(psiN_)

  end function geqdsk_equi2d_FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function geqdsk_equi2d_FdF(this, r) result(FdF)
  !
  ! return poloidal current [m T] and its first derivative
  !
  class(geqdsk_equi2d), intent(in) :: this
  real(real64),         intent(in) :: r(:)
  real(real64)                     :: FdF(0:1)

  real(real64) :: psiN


  psiN = (this%Psi%eval(r(1:2)) - this%Simag) / (this%Sibry - this%Simag)
  if (psiN > 1.d0  .or.  psiN < 0.d0) psiN = 1.d0
  FdF = this%Ffunc%deriv(psiN, 1)
  if (psiN == 1.d0) then
     FdF(1) = 0.d0
  else
     FdF(1) = FdF(1) / (this%Sibry - this%Simag)
  endif

  end function geqdsk_equi2d_FdF
  !-----------------------------------------------------------------------------
! type geqgsk_equi2d ===========================================================



! type sonnet_equi2d ===========================================================
! constructors:
  !-----------------------------------------------------------------------------
  function load_sonnet_equi2d(filename, scale_Bt, scale_Ip) result(this)
  !
  ! load equilibrium in sonnet format
  !
  use ieee_arithmetic
  use moose_error
  use moose_utils,    only: dirname
  use moose_analysis, only: bspline2d
  use flare_bfield,   only: init_magnetic_field
  character(len=*),     intent(in) :: filename
  real(real64),         intent(in), optional :: scale_Bt, scale_Ip
  type(sonnet_equi2d)              :: this

  real(real64), allocatable :: R(:), Z(:), Psi(:,:)
  character(len=120)        :: str
  real(real64) :: btf, rtf, psib, psix
  integer      :: i, iu, jm, km


  open  (newunit=iu, file=filename, action="read")
  ! skip past header lines
  do i=1,10
     read (iu, 1000) str
  enddo

  ! read grid resolution
  read  (iu, 1000) str;  read (str(18:26), *) jm
  read  (iu, 1000) str;  read (str(18:26), *) km
  print 3001, jm, km


  ! allocate memory for data
  allocate (R(jm), Z(km), Psi(jm,km))


  ! read characteristic parameters
  read  (iu, 1000) str;  read (str(14:40), *) psib
  read  (iu, 1000) str;  read (str(14:40), *) btf
  read  (iu, 1000) str;  read (str(14:40), *) rtf
  if (present(scale_Bt)) btf = btf * scale_Bt
  print 3003, btf
  print 3004, rtf


  ! read grid nodes [m]
  do i=1,2;  read (iu, 1000) str;  enddo
  read  (iu, *) R
  do i=1,2;  read (iu, 1000) str;  enddo
  read  (iu, *) Z


  ! read magnetic flux at grid points [Wb/rad]
  do i=1,2;  read (iu, 1000) str;  enddo
  read  (iu, *) Psi;  Psi = Psi + psib
  close (iu)
  if (present(scale_Ip)) Psi = Psi * scale_Ip
 1000 format(a120)
 3001 format(8x,'Grid resolution: ',19x,'R : ',i4,' points;   Z : ',i4,' points')
 3003 format(8x,'Vacuum toroidal magnetic field at R0:      ',2x,f8.3," T")
 3004 format(8x,'Reference position R0:                     ',2x,f8.3," m")


  ! sanity check
  if (any(ieee_is_nan(R))) call ERROR("invalid R", "load_sonnet_equi2d")
  if (any(ieee_is_nan(Z))) call ERROR("invalid Z", "load_sonnet_equi2d")
  if (any(ieee_is_nan(Psi))) call ERROR("invalid Psi", "load_sonnet_equi2d")

  ! initialize equi2d
  call init_magnetic_field(this, R(1), R(jm), Z(1), Z(km))
  this%F0 = btf * rtf
  allocate (this%Psi, source=bspline2d(R, Z, Psi))

  this%delta_Psi = psib - this%Psi%eval([rtf, 0.d0])
  call setup(this, rtf, 0.d0, dirname(filename))

  end function load_sonnet_equi2d
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine sonnet_equi2d_broadcast(this)
  use moose_mpi
  use moose_analysis, only: bspline2d
  class(sonnet_equi2d), intent(inout) :: this


  if (rank > 0) allocate (bspline2d :: this%Psi)
  call this%equi2d_broadcast()
  call proc(0)%broadcast(this%F0)

  end subroutine sonnet_equi2d_broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function sonnet_equi2d_FpsiN(this, psiN) result(F)
  !
  ! return toroidal field function [T m]
  !
  class(sonnet_equi2d), intent(in) :: this
  real(real64),         intent(in) :: psiN
  real(real64)                     :: F


  F = this%F0

  end function sonnet_equi2d_FpsiN
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function sonnet_equi2d_FdF(this, r) result(FdF)
  !
  ! return toroidal field function [T m] and its first derivative with respect to psiN
  !
  class(sonnet_equi2d), intent(in) :: this
  real(real64),         intent(in) :: r(:)
  real(real64)                     :: FdF(0:1)


  FdF(0) = this%F0
  FdF(1) = 0.d0

  end function sonnet_equi2d_FdF
  !-----------------------------------------------------------------------------
! type sonnet_equi2d ===========================================================

end module flare_equi2d
