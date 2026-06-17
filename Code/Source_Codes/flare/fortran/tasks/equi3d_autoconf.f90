subroutine flare_task_equi3d_autoconf(R0, Z0, Phi0, nsym, nphi, refinement)
  !
  ! Automatic configuriton of non-symmetric equilibrium.
  !
  ! **Parameters:**
  !
  ! :R0:    Initial guess of major radius [m] of magnetic axis at Phi0.
  !
  ! :Z0:    Initial guess of vertical position [m] of magnetic axis at Phi0.
  !
  ! :Phi0:  Toroidal position [deg] for initial guess of magnetic axis.
  !
  ! :nsym:  Toroidal symmetry of magnetic axis (default: symmetry of magnetic field).
  !
  ! :nphi:  Toroidal resolution of magnetic axis (default: 360 / nsym).
  !
  ! :refinement:  Automatic refinement of magnetic axis position.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_math, only: pi
  use flare_model
  use flare_poincare_map
  use flare_control
  use flare_tasks
  implicit none
  real(real64), intent(in) :: R0, Z0, Phi0
  integer,      intent(in) :: nsym, nphi
  logical,      intent(in) :: refinement

  integer, parameter :: nr = 31, nz = 31, ntransits = 10, nrefine_max = 10

  real(real64), parameter :: epsabs = 1.d-3

  type(poincare_map), allocatable :: P(:)
  real(real64) :: R00, Z00, Phi00, dr, dz, err
  character(len=128) :: filename
  integer :: i, iu, n, m


  R00   = R0
  Z00   = Z0
  Phi00 = Phi0/180.d0*pi
  ! guess R00 from equilibrium domain
  if (R00 <= 0.d0) R00 = (bfield%equi%lb(1) + bfield%equi%ub(1)) / 2


  ! display task parameters
  m = nsym;   if (nsym == 0) m = bfield%nfp
  call begin_task()
  if (rank == 0) then
     print *, "Automatic configuration of non-symmetric equilibrium ..."
     print *
     if (refinement) then
        print 1001, R00, Z00, Phi0
     else
        print 1002, R00, Z00, Phi0
     endif
     print *
     print 1003, m
     print *
  endif
 1001 format(3x,"- Initial guess for magnetic axis: (",f0.4," m, ",f0.4," m, ",f0.4," deg)")
 1002 format(3x,"- Reference point on magnetic axis: (",f0.4," m, ",f0.4," m, ",f0.4," deg)")
 1003 format(3x,"- Toroidal symmetry: ",i0)


  ! automatic refinement of magnetic axis position
  if (refinement) then
     if (rank == 0) print 2000
     dr = R00 / 20;   dz = dr
     do i=1,nrefine_max
        call search_grid(R00, Z00, dr, dz, err)
        if (rank == 0) print 2001, i, R00, Z00, err
        if (err < epsabs) exit
        if (i == nrefine_max  .and.  rank == 0) print *, "reached max. number of refinement steps"
     enddo
     if (rank == 0) print *
  endif
 2000 format(3x,"- Iterative approximation of magnetic axis position")
 2001 format(8x,"step ",i0,": magnetic axis = (",f0.6,", ",f0.6,"), accuracy = ",f0.6)


  ! construct magnetic axis from average of Poincare map
  n = nphi;   if (n <= 0) n = 360 / m
  if (rank == 0) then
     if (rank == 0) print 3000
     allocate (P(0:n-1), source=poincare_maps([R00, Z00, Phi00], 1, 0.d0, m, 1024, n))
     ! output Poincare maps
     if (verbose) then
     do i=0,n-1
        write (filename, 3002) i
        open  (newunit=iu, file=filename)
        write (iu, '(dt)') P(i)
        close (iu)
     enddo
     endif

     ! output magnetic axis
     open  (newunit=iu, file=".equi3d")
     write (iu, *) n, m
     do i=0,n-1
        R00 = sum(P(i)%points%column(1)) / P(i)%points%nelements()
        Z00 = sum(P(i)%points%column(2)) / P(i)%points%nelements()
        write (iu, *) R00, Z00, 360.d0 / m / n * i
     enddo
     close (iu)
  endif
  call finalize_task()
 3000 format(3x,"- Constructing magnetic axis position from Poincare map")
 3002 format("magnetic_axis_",i0,".dat")

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine search_grid(R00, Z00, dr, dz, err)
  use moose_error
  use moose_math, only: pi2
  use moose_grids, only: rmesh
  use flare_fieldline, only: fdriver
  use flare_control
  real(real64), intent(inout) :: R00, Z00
  real(real64), intent(inout) :: dr, dz
  real(real64), intent(  out) :: err

  real(real64), allocatable :: dist(:)
  type(fdriver) :: F
  type(rmesh)   :: mesh
  real(real64)  :: y(3), rz0(2)
  integer       :: i, imin(1), istat, j, k(2), n


  mesh = rmesh(R00-dr/2, R00+dr/2, nr, Z00-dz/2, Z00+dz/2, nz)
  n = mesh%nnodes()
  F = fdriver()
  call progress_bar(0, n)
  allocate (dist(0:n-1), source=0.d0)
  do i=rank,n-1,nproc
     y(1:2) = mesh%node(i)
     y(3)   = Phi00
     call F%reset()

     do j=1,ntransits
        if (very_verbose) print *, j, y
        rz0   = y(1:2)
        istat = F%evolve(y(3), y(3) + pi2/m, y(1:2))
        if (istat == SUCCESS) then
           dist(i) = dist(i) + sqrt(sum((y(1:2)-rz0)**2))
        else
           dist(i) = huge(1.d0)
           exit
        endif
     enddo
     call progress_bar(i+1, n)
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(dist)

  ! update position of magnetic axis and next search intervals
  imin = minloc(dist)-1
  k    = mesh%node_index(imin(1))
  rz0  = mesh%node(imin(1))
  dr   = max(2 * abs(rz0(1) - R00), epsabs);   if (k(1) == 0  .or.  k(1) == nr-1) dr = dr * 2
  dz   = max(2 * abs(rz0(2) - Z00), epsabs);   if (k(2) == 0  .or.  k(2) == nz-1) dz = dz * 2
  R00  = rz0(1);   Z00 = rz0(2)
  err  = dist(imin(1))

  call F%free()
  deallocate (dist)

  end subroutine search_grid
  !-----------------------------------------------------------------------------

end subroutine flare_task_equi3d_autoconf
