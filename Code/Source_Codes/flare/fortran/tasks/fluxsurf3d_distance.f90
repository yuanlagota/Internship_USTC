subroutine flare_task_fluxsurf3d_distance(fluxsurf3d_grid, nr, nz, rmin, rmax, zmin, zmax, output)
  !
  ! Compute R-Z distance to flux surface on cylindrical mesh.
  !
  ! **Parameters:**
  !
  ! :fluxsurf3d_grid:  Filename for flux surface discretization.
  !
  ! :nr, nz:           Number of points in R- and Z-direction.
  !
  ! :rmin, rmax:       Lower and upper limits in R-direction [m].
  !
  ! :zmin, zmax:       Lower and upper limits in Z-direction [m].
  !
  ! :output:           Filename for output of 3D B-Spline.
  !
  use iso_fortran_env
  use moose_math
  use moose_grids
  use moose_analysis, only: bspline3d
  use moose_geometry, only: polygon2d
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: fluxsurf3d_grid, output
  integer,          intent(in) :: nr, nz
  real(real64),     intent(in) :: rmin, rmax, zmin, zmax

  type(tpzmesh3d) :: M
  type(polygon2d) :: P
  type(bspline3d) :: B
  real(real64), allocatable :: r(:), z(:), phi(:), d(:,:,:)
  integer :: i, j, k, l, n, nphi


  call begin_task()
  if (report) then
     print *, "Computing R-Z distance to flux surface on cylindrical mesh:"
     print *

     M = tpzmesh3d(fluxsurf3d_grid)
     print 1000, trim(fluxsurf3d_grid)
     print 1001, M%n(1)
     print *

     print 1010, nr, nz
     print 1011, rmin, rmax
     print 1012, zmin, zmax
     print *
  endif
  call M%broadcast()
  nphi = M%n(1)
 1000 format(3x,"- Flux surface grid: ",a)
 1001 format(8x,"Number of toroidal points: ",i0)
 1010 format(3x,"- R-Z mesh with ",i0," x ",i0," points:")
 1011 format(8x,"R-range: ",f8.3," m -> ",f8.3," m")
 1012 format(8x,"Z-range: ",f8.3," m -> ",f8.3," m")


  ! construct domain for sampling
  allocate (r, source=linspace(rmin, rmax, nr))
  allocate (z, source=linspace(zmin, zmax, nz))
  allocate (phi, source=M%domain%u / 180.d0 * pi)


  ! sample R-Z distance on mesh
  allocate (d(nr, nz, nphi), source=0.d0)
  n = nr * nz * nphi
  l = 0
  call progress_bar(0, n)
  do k=1,nphi
     P = polygon2d(M%x(:,k-1,:))
     do i=1,nr
     do j=1,nz
        l = l + 1
        if (mod(l,nproc) /= rank) cycle
        d(i,j,k) = P%get_distance([r(i), z(j)])
        call progress_bar(l, n)
     enddo
     enddo
     call P%free()
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(d)


  ! compute B-spline coefficients
  if (rank == 0) then
     B = bspline3d(r, z, phi, d)
     call B%savenc(output)
  endif


  call finalize_task()

end subroutine flare_task_fluxsurf3d_distance
