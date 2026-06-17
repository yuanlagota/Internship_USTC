subroutine flare_task_magnetic_field(grid, dformat, output)
  !
  ! Sample magnetic field on a grid.
  !
  ! **Parameters:**
  !
  ! :grid:    Filename for grid with sample positions.
  !
  ! :dformat: Data format:
  !
  !           :-1: BR, BZ, Bphi (total field)
  !
  !           :0:  BR0, BZ0, Bphi0 (equilibrium field only)
  !
  !           :1:  BR1, BZ1, Bphi1 (perturbation field only)
  !
  !           :2:  Bmod, BR, BZ, Bphi, psiN
  !
  !           :3:  Bpsi1, Bpol1, Bphi1
  !
  !           :4:  psi, dpsi/dr, dpsi/dz, H11, H12, H22
  !
  ! :output:  Filename for output of data set.
  !
  use iso_fortran_env
  use moose_mpi
  use moose_dataset
  use moose_grids,  only: r3grid
  use flare_control
  use flare_model,  only: bfield
  use flare_equi2d, only: equi2d
  use flare_tasks
  implicit none
  character(len=*), intent(in) :: grid, output
  integer,          intent(in) :: dformat

  real(real64), pointer :: element(:)
  type(r3grid)  :: G
  type(dataset) :: D
  real(real64)  :: x(3), B(3), Bmod, B0(3), B1(3), Bpsi1, Bpol, Bpol1, psiN
  real(real64)  :: psi, dpsi(2), H(2,2)
  integer       :: i, n


  call begin_task()
  if (rank == 0) then
     print *, "Sampling magnetic field on grid ", trim(grid)
     print *

     G = r3grid(grid)
     print 1000, G%nnodes()
     print *
  endif
  call G%broadcast()
 1000 format(3x,"- Grid resolution: ",i0," nodes")


  ! set output format
  n = G%nnodes()
  select case(dformat)
  ! -1: BR, BZ, Bphi
  case(-1)
     D = dataset(3, G%nnodes())
     call D%set_metadata(1, 'BR',    "BR",   "T")
     call D%set_metadata(2, 'BZ',    "BZ",   "T")
     call D%set_metadata(3, 'Bphi',  "Bphi", "T")
     call D%set_expression("Bmod", "sqrt(BR**2 + BZ**2 + Bphi**2)", units="T")

  ! 0: BR, BZ, Bphi (equilibrium field)
  case(0)
     D = dataset(3, G%nnodes())
     call D%set_metadata(1, 'BR0',    "BR0",   "T")
     call D%set_metadata(2, 'BZ0',    "BZ0",   "T")
     call D%set_metadata(3, 'Bphi0',  "Bphi0", "T")
     call D%set_expression("Bmod", "sqrt(BR0**2 + BZ0**2 + Bphi0**2)", units="T")

  ! 1: BR, BZ, Bphi (perturbation field)
  case(1)
     D = dataset(3, G%nnodes())
     call D%set_metadata(1, 'BR1',   "BR1",   "T")
     call D%set_metadata(2, 'BZ1',   "BZ1",   "T")
     call D%set_metadata(3, 'Bphi1', "Bphi1", "T")
     call D%set_expression("Bmod", "sqrt(BR1**2 + BZ1**2 + Bphi1**2)", units="T")

  ! 2: |B|, BR, BZ, Bphi, psiN
  case(2)
     D = dataset(5, n)
     call D%set_metadata(1, 'Bmod',  "Bmod", "T")
     call D%set_metadata(2, 'BR',    "BR",   "T")
     call D%set_metadata(3, 'BZ',    "BZ",   "T")
     call D%set_metadata(4, 'Bphi',  "Bphi", "T")
     call D%set_metadata(5, 'psiN',  "Normalized Poloidal Flux")

  ! 3: Bpsi1, Bpol1, Bphi1
  case(3)
     D = dataset(3, n)
     call D%set_metadata(1, 'Bpsi1',  "Bpsi1", "T")
     call D%set_metadata(2, 'Bpol1',  "Bpol1", "T")
     call D%set_metadata(3, 'Bphi1',  "Bphi1", "T")

  ! 4:
  case(4)
     D = dataset(6, n)
     call D%set_metadata(1, 'psi',    "psi",    "Wb")
     call D%set_metadata(2, 'dpsidr', "dpsidr", "Wb/m")
     call D%set_metadata(3, 'dpsidz', "dpsidz", "Wb/m")
     call D%set_metadata(4, 'H11',    "H11",    "Wb/m**2")
     call D%set_metadata(5, 'H12',    "H12",    "Wb/m**2")
     call D%set_metadata(6, 'H22',    "H22",    "Wb/m**2")

  case default
     print 9000, dformat;   stop
  end select
 9000 format("ERROR: invalid output format ",i0)


  ! main loop ..........................................................
  call progress_bar(0, n)
  do i=rank,n-1,nproc
     x  = G%node(i)
     if (verbose) print *, i, x
     if (bfield%out_of_bounds(x)) cycle

     ! evaluate magnetic field
     B0   = bfield%equi%eval(x)
     B1   = bfield%perturbation_eval(x)
     B    = B0 + B1
     Bmod = sqrt(sum(B**2))
     Bpol = sqrt(sum(B(1:2)**2))
     psiN = bfield%equi%psiN(x)


     ! store results
     element => D%element(i)
     select case(dformat)
     case(-1)
        element = B

     case(0)
        element = B0

     case(1)
        element = B1

     case(2)
        element = [Bmod, B, psiN]

     case(3)
        Bpsi1 = - (B1(1)*B0(2) - B1(2)*B0(1)) / Bpol * bfield%equi%Bp_sign
        Bpol1 = sum(B1(1:2)*B0(1:2)) / Bpol
        element = [Bpsi1, Bpol1, B1(3)]

     case(4)
        select type(E => bfield%equi)
        class is (equi2d)
           psi = E%Psi%eval(x(1:2))
           dpsi = E%Psi%deriv(x(1:2))
           H = E%Psi%hessian(x(1:2))
           element = [psi, dpsi, H(1,1), H(1,2), H(2,2)]
        end select

     end select
     call progress_bar(i+1, n)
  enddo
  call finalize_progress_bar()
  call D%allreduce()


  if (rank == 0) then
     print 2000, trim(output)
     call D%set_geometry(grid, output)
     call D%savetxt(output)
  endif
  call D%free()
 2000 format(3x,"- Saving results to: ",a)
  call finalize_task()

end subroutine flare_task_magnetic_field
