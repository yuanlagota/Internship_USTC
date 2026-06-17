module flare_mmesh_utils
  use iso_fortran_env
  use flare_fieldline, only: fdriver
  implicit none
  private


  ! workspace for construction of field line segments
  type, extends(fdriver), public :: mdriver
     logical :: report

     contains
     procedure :: trace_vertices
  end type mdriver



  public :: &
     magnetic_flux_cells, magnetic_flux_xsects

  contains
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine trace_vertices(this, n, x0, phi, it0, x, iend)
  !
  ! Trace vertices across toroidal domain
  !
  ! Input:
  !    n           number of vertices
  !    x0          initial (r,z) coordinates
  !    phi         toroidal positions [rad]
  !    it0         index for initial toroidal position in phi-array
  !
  ! Output:
  !    x(1:2),*)   (r,z) coordinates of vertices at phi locations
  !    x(3:4),*)   (br/bphi, bz/bphi) values at x(1:2,*)
  !    x(5,*)      magnetic field strength at x(1:2,*)
  !    iend(1,*)   final index for field line tracing in backward direction
  !    iend(2,*)   final index for field line tracing in forward direction
  !
  use moose_mpi
  use flare_control, only: progress_bar, finalize_progress_bar
  use flare_model,   only: bfield
  class(mdriver), intent(inout) :: this
  integer,        intent(in   ) :: n, it0
  real(real64),   intent(in   ) :: x0(2, *), phi(0:)
  real(real64),   intent(  out) :: x(5, 0:size(phi)-1, *)
  integer,        intent(  out) :: iend(2, *)

  integer, parameter :: index_for_kdir(-1:1) = [1, 0, 2]

  real(real64) :: b(3), phi0, y(3)
  integer :: i, istat, iu, k, kdir, kend, nphi


  nphi = size(phi)
  phi0 = phi(it0)
  x(:, 0:nphi-1, 1:n) = 0.d0
  iend(:, 1:n) = it0


  if (this%report) call progress_bar(0, n)
  do i=1+rank,n,nproc
     ! initial position
     b = bfield%eval([x0(:,i), phi0])
     x(:, it0, i) = [x0(:,i), b(1:2) / b(3), norm2(b)]

     ! trace in forward and backward direction across domain
     do kdir=-1,1,2
        call this%reset()
        y = [x0(:,i), phi0]

        kend = 0;   if (kdir == 1) kend = nphi-1
        do k=it0+kdir,kend,kdir
           istat = this%evolve3(y, phi(k))
           if (istat > 0) exit

           b = bfield%eval(y)
           x(:, k, i) = [y(1:2), b(1:2) / b(3), norm2(b)]
           iend(index_for_kdir(kdir), i) = k
        enddo
     enddo
     if (this%report) call progress_bar(i, n)
  enddo
  if (this%report) call finalize_progress_bar()


  call moose_mpi_sum(x(:,0:nphi-1,1:n))
  call moose_mpi_sum(iend(:,1:n))

  end subroutine trace_vertices
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  pure function magnetic_flux_cells(phi, f) result(flux)
  !
  ! compute magnetic flux in cells of flux tube
  !
  use moose_uqwork, only: area
  real(real64), intent(in) :: phi(:)
  real(real64), intent(in) :: f(5, size(phi), 4)
  real(real64)             :: flux(2:size(phi))

  real(real64) :: a1, a2, b1, b2, c1(2), c2(2), dfl, pitch
  integer :: k, n


  n = size(phi)
  do k=1,n
     a2 = area(f(1:2,k,1), f(1:2,k,2), f(1:2,k,3), f(1:2,k,4))
     b2 = sum(f(5,k,:)) / 4
     c2 = sum(f(1:2,k,:), dim=2) / 4

     if (k > 1) then
        dfl = (c1(1) + c2(1)) / 2 * abs(phi(k) - phi(k-1))
        pitch = dfl / sqrt(dfl**2 + sum((c2 - c1)**2))
        flux(k) = (a1 + a2) * pitch * (b1 + b2) / 4
     endif

     a1 = a2
     b1 = b2
     c1 = c2
  enddo

  end function magnetic_flux_cells
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function magnetic_flux_xsects(phi, f) result(flux)
  !
  ! compute magnetic flux at cross-sections of flux tube
  !
  use moose_uqwork, only: area
  real(real64), intent(in) :: phi(:)
  real(real64), intent(in) :: f(5, size(phi), 4)
  real(real64)             :: flux(size(phi))

  real(real64) :: a, b, g(2), pitch
  integer :: k, n


  n = size(phi)
  do k=1,n
     a = area(f(1:2,k,1), f(1:2,k,2), f(1:2,k,3), f(1:2,k,4))
     b = sum(f(5,k,:)) / 4
     g = sum(f(3:4,k,:), dim=2) / 4

     pitch = 1.d0 / sqrt(1.d0  +  g(1)**2  +  g(2)**2)
     flux(k) = a * b * pitch
  enddo

  end function magnetic_flux_xsects
  !-----------------------------------------------------------------------------

end module flare_mmesh_utils
