subroutine flare_task_grazing_angle(surface_mesh, output)
  !
  ! Compute field line grazing angle on surface mesh
  !
  ! **Parameters:**
  !
  ! :surface_mesh:   Filename for surface mesh (r3grid with tpzmesh3d domain).
  !
  ! :output:         Filename for output.
  !
  use iso_fortran_env
  use moose_error
  use moose_mpi
  use moose_math,     only: pi
  use moose_grids,    only: r3grid, tpzmesh3d
  use moose_geometry, only: torosurf
  use moose_data,     only: dataset
  use flare_control,  only: progress_bar, finalize_progress_bar
  use flare_model
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: surface_mesh, output

  type(r3grid)  :: G
  type(tpzmesh3d), pointer :: M
  type(torosurf) :: T
  type(dataset) :: D
  real(real64) :: b(3), x(3), v(3)
  integer :: i, k(2), n


  call begin_task()
  if (rank == 0) then
     print *, "Computing field line grazing angle on surface mesh ", trim(surface_mesh)
     print *

     G = r3grid(surface_mesh)
     if (any(G%map /= [1, 2, 3])) call ERROR("invalid map for r3grid domain")
  endif
  call G%broadcast()


  ! verify tpzmesh3d domain and convert units to m, rad
  select type(domain => G%domain)
  class is (tpzmesh3d)
     M => domain
     M%x = M%x * G%length_scale
     if (G%in_degrees) M%domain%u = M%domain%u * pi / 180

  class default
     call ERROR("invalid r3grid domain type")
  end select
  T = torosurf(M)


  ! main loop
  n = M%ncells()
  D = dataset(1, n)
  call D%set_metadata(1, 'alpha',  "Grazing angle", "deg")
  call progress_bar(0, n)
  do i=rank,n-1,nproc
     k = M%cell_index(i)

     ! evaluate magnetic field vector and surface normal at cell center
     x = T%interp(k + 0.5d0)
     v = T%normal_vector(k + 0.5d0)
     b = bfield%eval(x);   b = b / norm2(b)
     D%values(1, i) = asin(sum(b * v)) / pi * 180

     call progress_bar(i+1, n)
  enddo
  call finalize_progress_bar()
  call D%allreduce()


  ! save results
  if (rank == 0) then
     call D%set_geometry(surface_mesh, output)
     call D%savetxt(output)
  endif
  call finalize_task()

end subroutine flare_task_grazing_angle
