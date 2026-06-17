module flare_interp
  use iso_fortran_env
  use flare_bfield
  implicit none
  private


  character(len=*), parameter :: &
     ASCII  = "ascii", &
     BINARY = "binary"


  type, extends(magnetic_field), public :: interp_bfield
     real(real64), allocatable :: b(:,:)
     integer      :: nr, nz, nphi
     real(real64) :: dr, dz, dphi, fp

     ! interpolation accelerators
     real(real64), pointer :: bo(:,:)
     integer, pointer :: li, lj, lk

     contains
     procedure :: broadcast
     procedure :: free

     procedure :: update_bo
     procedure :: eval
     procedure :: jac

     procedure :: export
  end type interp_bfield


  interface interp_bfield
     procedure :: sample_bfield
  end interface interp_bfield



  public :: &
     new_interp_bfield, &
     load_interp_bfield

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  subroutine aux_init_interp_bfield(this)
  use moose_math, only: pi2
  class(interp_bfield), intent(inout) :: this

  integer :: m


  m       = this%nr * this%nz * this%nphi
  this%fp = pi2 / this%nfp
  allocate (this%b(8, 0:m-1), source=0.d0)
  allocate (this%li, source=-1)
  allocate (this%lj, source=-1)
  allocate (this%lk, source=-1)
  allocate (this%bo(20,3), source=0.d0)

  end subroutine aux_init_interp_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function new_interp_bfield(nr, nz, nphi, nfp, rmin, rmax, zmin, zmax) result(this)
  integer,      intent(in) :: nr, nz, nphi, nfp
  real(real64), intent(in) :: rmin, rmax, zmin, zmax
  type(interp_bfield)      :: this


  this%nr = nr
  this%nz = nz
  this%nphi = nphi
  call init_magnetic_field(this, rmin, rmax, zmin, zmax, nfp)
  call aux_init_interp_bfield(this)
  this%dr = (rmax - rmin) / (nr - 1)
  this%dz = (zmax - zmin) / (nz - 1)
  this%dphi = this%fp / nphi

  end function new_interp_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function sample_bfield(bfield, nr, nz, nphi, nfp, rmin, rmax, zmin, zmax) result(this)
  use moose_math
  use moose_mpi
  use flare_control
  class(magnetic_field), intent(in) :: bfield
  integer,               intent(in) :: nr, nz, nphi, nfp
  real(real64),          intent(in) :: rmin, rmax, zmin, zmax
  type(interp_bfield)               :: this

  real(real64) :: b(3), jac(3,3), r(nr), z(nz), phi(nphi), x(3)
  integer :: i, j, k, l, n


  this = new_interp_bfield(nr, nz, nphi, nfp, rmin, rmax, zmin, zmax)
  n = nr * nz * nphi
  r = linspace(rmin, rmax, nr)
  z = linspace(zmin, zmax, nz)
  phi = linspace(0.d0, pi2 / nfp, nphi)
  if (report) call progress_bar(0, n)
  do j=1,nz
  do i=1,nr
  do k=1,nphi
     l = (k-1) + ((i-1) + (j-1)*nr)*nphi
     if (.not. mod(l,nproc) == rank) cycle

     x = [r(i), z(j), phi(k)]
     b = bfield%eval(x)
     jac = bfield%jac(x)
     this%b(1,l) = b(3)
     this%b(2,l) = b(1)
     this%b(3,l) = b(2)
     this%b(4,l) = jac(3,1)
     this%b(5,l) = jac(3,2)
     this%b(6,l) = jac(1,1)
     this%b(7,l) = jac(1,2)
     this%b(8,l) = jac(2,2)
     if (report) call progress_bar(l+1, n)
  enddo
  enddo
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(this%b)

  end function sample_bfield
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load_interp_bfield(filename, filetype, amplitude, bfield_scale, length_scale) result(this)
  use moose_utils,  only: basename
  use flare_bfield, only: init_magnetic_field
  character(len=*), intent(in) :: filename, filetype
  real(real64),     intent(in) :: amplitude
  real(real64),     intent(in), optional :: bfield_scale, length_scale
  type(interp_bfield)          :: this

  character(len=72) :: title
  real(real32), allocatable :: b32(:,:)
  real(real64) :: domain(5), rdummy
  integer :: iu, idummy, j, nfp, nhellu, npollu, npaqlu, nspilu


  print *
  print 1000, basename(filename)
 1000 format(3x,"- Interpolation of magnetic field on cylindrical grid: ",a)

  ! open data file
  select case(filetype)
  case(ASCII)
     open (newunit=iu, file=filename, action='read')
  case(BINARY)
     open (newunit=iu, file=filename, action='read', form='unformatted')
  case default
     print 9000, trim(filetype);   stop
  end select
 9000 format("ERROR: unkown file type ",a)


  ! read header
  select case(filetype)
  case(ASCII)
     read (iu, '(a)') title
     read (iu, *) this%dphi, this%dr, this%dz, domain, this%nphi, this%nr, nfp, this%nz, &
                  nhellu, npollu, npaqlu, nspilu
  case(BINARY)
     read (iu) title, this%dphi, this%dr, this%dz, domain, this%nphi, this%nr, nfp, this%nz, &
               nhellu, npollu, npaqlu, nspilu
  end select
  if (present(length_scale)) then
     this%dr = this%dr * length_scale
     this%dz = this%dz * length_scale
     domain  = domain  * length_scale
  endif
  call init_magnetic_field(this, domain(1), domain(2), domain(4), domain(5), nfp)
  call aux_init_interp_bfield(this)
  print 1001, trim(adjustl(title))
  print 1002, amplitude
  print 1003, domain(1:2)
  print 1004, domain(4:5)
  print 1005, nfp
  print 1006, this%nphi, this%nr, this%nz
 1001 format(8x,a)
 1002 format(8x,"amplitude = ",f0.3)
 1003 format(8x,"R range: ",f0.3," -> ",f0.3," m")
 1004 format(8x,"Z range: ",f0.3," -> ",f0.3," m")
 1005 format(8x,"Symmetry: ",i0)
 1006 format(8x,"Resolution (nphi x nr x nz): ",i0," x ",i0," x ",i0)


  ! legacy input .......................................................
  if (nhellu > 0) then
     select case(filetype)
     case(ASCII)
        read (iu, *) (rdummy, j=1,nhellu), (rdummy, j=1,10), (idummy, j=1,nhellu+3)
     case(BINARY)
        read (iu) (rdummy, j=1,nhellu), (rdummy, j=1,10), (idummy, j=1,nhellu+3)
     end select
  endif
  if (npollu > 0) then
     select case(filetype)
     case(ASCII)
        read (iu, *) (idummy, j=1,20), (rdummy, j=1,20), (rdummy, j=1,21*20), &
                     (rdummy, j=1,21*20), (rdummy, j=1,21*20)
     case(BINARY)
        read (iu)    (idummy, j=1,20), (rdummy, j=1,20), (rdummy, j=1,21*20), &
                     (rdummy, j=1,21*20), (rdummy, j=1,21*20)
     end select
  endif
  if (npaqlu > 0) then
     select case(filetype)
     case(ASCII)
        read (iu, *) (idummy, j=1,10), (rdummy, j=1,10*10), (rdummy, j=1,10*10)
     case(BINARY)
        read (iu)    (idummy, j=1,10), (rdummy, j=1,10*10), (rdummy, j=1,10*10)
     end select
  endif
  if (nspilu > 0) then
     select case(filetype)
     case(ASCII)
        read (iu, *) (rdummy, j=1,7*nspilu)
     case(BINARY)
        read (iu)    (rdummy, j=1,7*nspilu)
     end select
  endif
  !.....................................................................


  ! read data
  select case(filetype)
  case(ASCII)
     read (iu, *) this%b
  case(BINARY)
     allocate (b32(8,0:ubound(this%b,2)))
     read (iu) b32;   this%b = b32
     deallocate (b32)
  end select
  close (iu)
  this%b = amplitude * this%b
  if (present(bfield_scale)) this%b = this%b * bfield_scale
  if (present(length_scale)) this%b(4:8,:) = this%b(4:8,:) / length_scale

  end function load_interp_bfield
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(interp_bfield), intent(inout) :: this


  call this%bfield_broadcast()
  call proc(0)%broadcast_allocatable(this%b)
  call proc(0)%broadcast(this%nr)
  call proc(0)%broadcast(this%nz)
  call proc(0)%broadcast(this%nphi)
  call proc(0)%broadcast(this%dr)
  call proc(0)%broadcast(this%dz)
  call proc(0)%broadcast(this%dphi)
  call proc(0)%broadcast(this%fp)
  if (rank > 0) then
     allocate (this%li, source=-1)
     allocate (this%lj, source=-1)
     allocate (this%lk, source=-1)
     allocate (this%bo(20,3), source=0.d0)
  endif

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(interp_bfield), intent(inout) :: this


  deallocate (this%b, this%li, this%lj, this%lk, this%bo)
  call this%mfunc_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine update_bo(this, i, j, k)
  class(interp_bfield), intent(in) :: this
  integer,              intent(in) :: i, j, k

  real(real64) :: bi(8,8), r1, r2
  integer      :: inode(8), k2, l


  if (k == this%lk  .and.  i == this%li  .and.  j == this%lj) return
  this%li  = i
  this%lj  = j
  this%lk  = k
  r1       = i * this%dr + this%lb(1)
  r2       = r1 + this%dr

  inode(1) = k  + (i    + j   *this%nr)*this%nphi
  inode(3) = k  + (i+1  + j   *this%nr)*this%nphi
  inode(5) = k  + (i    +(j+1)*this%nr)*this%nphi
  inode(7) = k  + (i+1  +(j+1)*this%nr)*this%nphi

  k2 = mod(k+1, this%nphi)
  inode(2) = k2 + (i    + j   *this%nr)*this%nphi
  inode(4) = k2 + (i+1  + j   *this%nr)*this%nphi
  inode(6) = k2 + (i    +(j+1)*this%nr)*this%nphi
  inode(8) = k2 + (i+1  +(j+1)*this%nr)*this%nphi

  do l=1,8
     bi(1:8,l) = this%b(l,inode)
  enddo

  this%bo(1:8,1:3)=  bi(1:8,1:3)
  this%bo( 9,1)=-r1*(bi(1,6)+bi(1,8)-bi(2,6)-bi(2,8))-bi(1,2)+bi(2,2)
  this%bo( 9,2)= r1*(bi(1,4)-bi(2,4))                +bi(1,1)-bi(2,1)
  this%bo( 9,3)= r1*(bi(1,5)-bi(2,5))
  this%bo(10,1)=     bi(1,4)-bi(3,4)
  this%bo(10,2)=     bi(1,6)-bi(3,6)
  this%bo(10,3)=     bi(1,7)-bi(3,7)
  this%bo(11,1)=     bi(1,5)-bi(5,5)
  this%bo(11,2)=     bi(1,7)-bi(5,7)
  this%bo(11,3)=     bi(1,8)-bi(5,8)
  this%bo(12,1)=-r2*(bi(3,6)+bi(3,8)-bi(4,6)-bi(4,8))-bi(3,2)+bi(4,2)
  this%bo(12,2)= r2*(bi(3,4)-bi(4,4))                +bi(3,1)-bi(4,1)
  this%bo(12,3)= r2*(bi(3,5)-bi(4,5))
  this%bo(13,1)=     bi(2,4)-bi(4,4)
  this%bo(13,2)=     bi(2,6)-bi(4,6)
  this%bo(13,3)=     bi(2,7)-bi(4,7)
  this%bo(14,1)=     bi(2,5)-bi(6,5)
  this%bo(14,2)=     bi(2,7)-bi(6,7)
  this%bo(14,3)=     bi(2,8)-bi(6,8)
  this%bo(15,1)=-r1*(bi(5,6)+bi(5,8)-bi(6,6)-bi(6,8))-bi(5,2)+bi(6,2)
  this%bo(15,2)= r1*(bi(5,4)-bi(6,4))                +bi(5,1)-bi(6,1)
  this%bo(15,3)= r1*(bi(5,5)-bi(6,5))
  this%bo(16,1)=     bi(5,4)-bi(7,4)
  this%bo(16,2)=     bi(5,6)-bi(7,6)
  this%bo(16,3)=     bi(5,7)-bi(7,7)
  this%bo(17,1)=     bi(3,5)-bi(7,5)
  this%bo(17,2)=     bi(3,7)-bi(7,7)
  this%bo(17,3)=     bi(3,8)-bi(7,8)
  this%bo(18,1)=-r2*(bi(7,6)+bi(7,8)-bi(8,6)-bi(8,8))-bi(7,2)+bi(8,2)
  this%bo(18,2)= r2*(bi(7,4)-bi(8,4))                +bi(7,1)-bi(8,1)
  this%bo(18,3)= r2*(bi(7,5)-bi(8,5))
  this%bo(19,1)=     bi(6,4)-bi(8,4)
  this%bo(19,2)=     bi(6,6)-bi(8,6)
  this%bo(19,3)=     bi(6,7)-bi(8,7)
  this%bo(20,1)=     bi(4,5)-bi(8,5)
  this%bo(20,2)=     bi(4,7)-bi(8,7)
  this%bo(20,3)=     bi(4,8)-bi(8,8)

  end subroutine update_bo
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function eval(this, x) result(B)
  class(interp_bfield), intent(in) :: this
  real(real64),         intent(in) :: x(this%ndim)
  real(real64)                     :: B(this%ndim)

  real(real64) :: r, z, phi, alr, alz, alp
  real(real64) :: rq(2), zq(2), pq(3), rq1, zq1, pq1, pr(4), pz(4), prz(20)
  integer      :: i, j, k


  r   = x(1);   z = x(2)
  phi = modulo(x(3), this%fp)
  if (r < this%lb(1)  .or.  r >= this%ub(1)  .or. &
      z < this%lb(2)  .or.  z >= this%ub(2)) then
     B = 0.d0
     return
  endif

  alr = (r  -this%lb(1)) / this%dr;     i = int(floor(alr))
  alz = (z  -this%lb(2)) / this%dz;     j = int(floor(alz))
  alp = (phi-this%lb(3)) / this%dphi;   k = int(floor(alp))
  rq(2) = alr - i;   rq(1) = 1 - rq(2)
  zq(2) = alz - j;   zq(1) = 1 - zq(2)
  pq(2) = alp - k;   pq(1) = 1 - pq(2)

  rq1   = rq(1) * rq(2) * this%dr/2
  zq1   = zq(1) * zq(2) * this%dz/2
  pq1   = pq(1) * pq(2) * this%dphi/2
  call this%update_bo(i, j, k)

  pr (1:2) = pq(1:2) * rq(1)
  pr (3:4) = pq(1:2) * rq(2)
  pz (1:2) = pq(1:2) * zq(1)
  pz (3:4) = pq(1:2) * zq(2)

  prz( 1:4) = pr(1:4) * zq (1)
  prz( 5:8) = pr(1:4) * zq (2)
  prz( 9  ) = rq(1) * zq(1) * pq1
  prz(12  ) = rq(2) * zq(1) * pq1
  prz(15  ) = rq(1) * zq(2) * pq1
  prz(18  ) = rq(2) * zq(2) * pq1
  prz(10  ) = pz(1) * rq1
  prz(13  ) = pz(2) * rq1
  prz(16  ) = pz(3) * rq1
  prz(19  ) = pz(4) * rq1
  prz(11  ) = pr(1) * zq1
  prz(14  ) = pr(2) * zq1
  prz(17  ) = pr(3) * zq1
  prz(20  ) = pr(4) * zq1

  B(3) = sum( prz(1:20)*this%bo(1:20,1) )
  B(1) = sum( prz(1:20)*this%bo(1:20,2) )
  B(2) = sum( prz(1:20)*this%bo(1:20,3) )

  end function eval
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function jac(this, x)
  class(interp_bfield), intent(in) :: this
  real(real64),         intent(in) :: x(this%ndim)
  real(real64)                     :: jac(this%mdim, this%ndim)

  real(real64) :: r, z, phi, alr, alz, alp
  real(real64) :: rq(2), zq(2), pq(3), rq1, zq1, pq1, pr(4), pz(4), prz(20), prz2(20), &
                  prz3(20), B(8), pq1dr1, pq1dz1, rq1dz1, zq1dr1, dr1, dz1, rq2, zq2
  integer      :: i, j, k


  r   = x(1);   z = x(2)
  phi = modulo(x(3), this%fp)

  alr = (r-this%lb(1)) / this%dr;     i = int(floor(alr))
  alz = (r-this%lb(2)) / this%dz;     j = int(floor(alz))
  alp = (r-this%lb(3)) / this%dphi;   k = int(floor(alp))
  dr1 = 1.d0 / this%dr;   dz1 = 1.d0 / this%dz
  rq(2) = alr - i;   rq(1) = 1 - rq(2)
  zq(2) = alz - j;   zq(1) = 1 - zq(2)
  pq(2) = alp - k;   pq(1) = 1 - pq(2)

  rq1   = rq(1) * rq(2) * this%dr/2
  zq1   = zq(1) * zq(2) * this%dz/2
  pq1   = pq(1) * pq(2) * this%dphi/2
  call this%update_bo(i, j, k)

  pr (1:2) = pq(1:2) * rq(1)
  pr (3:4) = pq(1:2) * rq(2)
  pz (1:2) = pq(1:2) * zq(1)
  pz (3:4) = pq(1:2) * zq(2)

  prz( 1:4) = pr(1:4) * zq (1)
  prz( 5:8) = pr(1:4) * zq (2)
  prz( 9  ) = rq(1) * zq(1) * pq1
  prz(12  ) = rq(2) * zq(1) * pq1
  prz(15  ) = rq(1) * zq(2) * pq1
  prz(18  ) = rq(2) * zq(2) * pq1
  prz(10  ) = pz(1) * rq1
  prz(13  ) = pz(2) * rq1
  prz(16  ) = pz(3) * rq1
  prz(19  ) = pz(4) * rq1
  prz(11  ) = pr(1) * zq1
  prz(14  ) = pr(2) * zq1
  prz(17  ) = pr(3) * zq1
  prz(20  ) = pr(4) * zq1

  B(3) = sum( prz(1:20)*this%bo(1:20,1) )
  B(1) = sum( prz(1:20)*this%bo(1:20,2) )
  B(2) = sum( prz(1:20)*this%bo(1:20,3) )

  prz2(1:2) = -dr1 * pz (1:2)
  prz2(3:4) = -prz2(1:2)
  prz2(5:6) = -dr1 * pz (3:4)
  prz2(7:8) = -prz2(5:6)

  prz3(1:2) = -dz1 * pr (1:2)
  prz3(5:6) = -prz3(1:2)
  prz3(3:4) = -dz1 * pr (3:4)
  prz3(7:8) = -prz3(3:4)

  pq1dr1   = -pq1*dr1
  prz2( 9) =  zq(1) * pq1dr1
  prz2(12) = -prz2( 9)
  prz2(15) =  zq(2) * pq1dr1
  prz2(18) = -prz2(15)
  pq1dz1   = -pq1*dz1
  prz3( 9) =  rq(1) * pq1dz1
  prz3(15) = -prz3( 9)
  prz3(12) =  rq(2) * pq1dz1
  prz3(18) = -prz3(12)
  rq1dz1   = -rq1*dz1
  prz3(10) =  pq(1) * rq1dz1
  prz3(16) = -prz3(10)
  prz3(13) =  pq(2) * rq1dz1
  prz3(19) = -prz3(13)
  zq1dr1   = -zq1*dr1
  prz2(11) =  pq(1) * zq1dr1
  prz2(17) = -prz2(11)
  prz2(14) =  pq(2) * zq1dr1
  prz2(20) = -prz2(14)
  rq2      = 0.5 - rq(2)
  prz2(10) =  rq2 * pz(1)
  prz2(13) =  rq2 * pz(2)
  prz2(16) =  rq2 * pz(3)
  prz2(19) =  rq2 * pz(4)
  zq2      = 0.5 - zq(2)
  prz3(11) =  zq2 * pr(1)
  prz3(14) =  zq2 * pr(2)
  prz3(17) =  zq2 * pr(3)
  prz3(20) =  zq2 * pr(4)

  B(4)     = sum( prz2(1:20)*this%bo(1:20,1) )
  B(5)     = sum( prz3(1:20)*this%bo(1:20,1) )
  B(6)     = sum( prz2(1:20)*this%bo(1:20,2) )
  B(7)     = sum( prz3(1:20)*this%bo(1:20,2) )
  B(8)     = sum( prz3(1:20)*this%bo(1:20,3) )

  ! dBr
  jac(1,1) = B(6)
  jac(1,2) = B(7)
  jac(1,3) = R*B(4) + B(3)
  ! dBz
  jac(2,1) = B(7)
  jac(2,2) = B(8)
  jac(2,3) = R*B(5)
  ! dBphi
  jac(3,1) = B(4)
  jac(3,2) = B(5)
  jac(3,3) = -R*(B(6)+B(8)) - B(1)
  jac      = jac

  end function jac
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine export(this, filename, filetype, title, bfield_scale, length_scale)
  use moose_error
  class(interp_bfield), intent(in) :: this
  character(len=*),     intent(in) :: filename, filetype, title
  real(real64),         intent(in), optional :: bfield_scale, length_scale

  character(len=72) :: s
  real(real64), allocatable :: b(:,:)
  real(real64) :: aux(7)
  integer :: iu, m


  aux(1) = this%dr
  aux(2) = this%dz
  aux(3) = this%lb(1)
  aux(4) = this%ub(1)
  aux(5) = (this%ub(1) + this%lb(1)) / 2
  aux(6) = this%lb(2)
  aux(7) = this%ub(2)

  m = this%nr * this%nz * this%nphi
  allocate (b(8, 0:m-1), source=this%b)

  ! user defined units
  if (present(bfield_scale)) b = b / bfield_scale
  if (present(length_scale)) then
     aux = aux / length_scale
     b(4:8,:) = b(4:8,:) * length_scale
  endif


  select case (filetype)
  case (ASCII)
     open  (newunit=iu, file=filename)
     write (iu, '(a72)') title
     write (iu, *) this%dphi, aux, this%nphi, this%nr, this%nfp, this%nz, 0, 0, 0, 0
     write (iu, *) b
     close (iu)

  case (BINARY)
     s = title
     open  (newunit=iu, file=filename, form='unformatted')
     write (iu) s
     write (iu) this%dphi, aux, this%nphi, this%nr, this%nfp, this%nz, 0, 0, 0, 0
     write (iu) real(b, real32)
     close (iu)

  case default
     call ERROR("invalid filetype = "//trim(filetype), "interp_bfield%export")
  end select
  deallocate (b)

  end subroutine export
  !-----------------------------------------------------------------------------

end module flare_interp
