#include <txtio.h>
!===============================================================================
! field line connection length analysis
!===============================================================================
module flare_connection_length
  use iso_fortran_env
  use moose_txtio
  implicit none
  private


  ! connection length histogram for field lines started on an equilibrium flux surface
  type, extends(txtio), public :: connection_histogram
     real(real64), allocatable :: psiN(:), turns(:)
     integer, allocatable :: counts(:,:,:)
     integer :: nsamples

     contains
     procedure :: free
     procedure :: cdf, pdf
     procedure :: write_formatted
  end type connection_histogram
  ! connection_histogram .......................................................



  public :: &
     compute_connection_histogram

  contains
  !-----------------------------------------------------------------------------


! auxiliary procedures:
  !-----------------------------------------------------------------------------
  subroutine aux_init_connection_histogram(this, psiN, turns)
  class(connection_histogram), intent(inout) :: this
  real(real64),                intent(in   ) :: psiN(:), turns(:)


  call init_txtio(this, "connection_histogram")
  allocate (this%psiN,  source=psiN)
  allocate (this%turns, source=turns)
  allocate (this%counts(0:size(turns), size(psiN), -1:1), source=0)

  end subroutine aux_init_connection_histogram
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine aux_scan(this, ipsiN, G)
  use moose_mpi
  use moose_math, only: deg3, pi2
  use moose_grids
  use flare_control
  use flare_fieldline
  class(connection_histogram), intent(inout) :: this
  integer,                     intent(in   ) :: ipsiN
  class(grid),                 intent(in   ) :: G

  type(fdriver) :: F
  real(real64)  :: x(3)
  integer       :: counts(0:size(this%turns), -1:1)
  integer       :: i, istat, j, jdir, n


  n = G%nnodes()
  F = fdriver()
  counts = 0.d0
  call progress_bar(0, n)
  do i=rank,n-1,nproc
     ! trace in both directions
     do jdir=-1,1,2
        call F%reset()

        x = G%node(i)
        do j=1,size(this%turns)
           istat = F%evolve3(x, jdir * pi2 * this%turns(j))
           if (istat == INTERSECT_BOUNDARY) then
              counts(j:, jdir) = counts(j:, jdir) + 1
              exit
           elseif (istat > 0) then
              print *, "initial point: ", deg3(G%node(i))
              call FIELDLINE_ERROR(F, istat)
           endif
        enddo
     enddo

     call progress_bar(i+1, n)
  enddo
  call finalize_progress_bar()
  call moose_mpi_sum(counts);   this%counts(:,ipsiN,:) = counts
  this%counts(:,:,0) = this%counts(:,:,-1) + this%counts(:,:,1)
  this%nsamples = n
  call F%free()

  end subroutine aux_scan
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function compute_connection_histogram(psiN, turns, nsym, nphi, ntheta, param) result(this)
  use moose_math
  use moose_grids
  use flare_control
  use flare_model
  use flare_fluxsurf2d
  real(real64),           intent(in) :: psiN(:), turns(:)
  integer,                intent(in) :: nsym, nphi, ntheta
  character(len=*),       intent(in) :: param
  type(connection_histogram)         :: this

  class(grid), allocatable :: G
  real(real64) :: x(2)
  type(cgrid) :: F
  integer     :: i


  call assert_equi2d("fluxsurf2d_connection_histogram")
  call aux_init_connection_histogram(this, psiN, turns)


  do i=1,size(psiN)
     if (report) print 1000, psiN(i)

     ! initial point for flux surface contour
     x = equi2d%rzcoords(psiN(i), 0.d0)

     ! construct sample points on flux surface
     F = fluxsurf2d_grid(equi2d%rzcoords(psiN(i), 0.d0), ntheta, endpoint=.false., param=param)
     allocate (G, source=rmesh3d(linspace(0.d0, pi2/nsym, nphi, endpoint=.false.), F))

     ! compute connection length histogram for i-th flux surface
     call aux_scan(this, i, G)
     if (report) print *

     ! cleanup
     call G%free();   deallocate (G)
  enddo
 1000 format(8x,"psiN = ",f0.4," flux surface")

  end function compute_connection_histogram
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(connection_histogram), intent(inout) :: this


  deallocate (this%psiN, this%turns, this%counts)
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function cdf(this, ipsiN, jdir)
  class(connection_histogram), intent(in) :: this
  integer,                     intent(in) :: ipsiN, jdir
  real(real64)                            :: cdf(0:size(this%turns))


  cdf = 1.d0 * this%counts(:,ipsiN,jdir) / this%nsamples

  end function cdf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function pdf(this, ipsiN, jdir)
  class(connection_histogram), intent(in) :: this
  integer,                     intent(in) :: ipsiN, jdir
  real(real64)                            :: pdf(size(this%turns))

  real(real64) :: cdf(0:size(this%turns)), turns(0:size(this%turns))
  integer :: n


  turns(0)  = 0.d0
  turns(1:) = this%turns

  n   = size(this%turns)
  cdf = this%cdf(ipsiN,jdir)
  pdf = (cdf(1:n) - cdf(0:n-1)) / (turns(1:) - turns(0:n-1))

  end function pdf
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(connection_histogram), intent(in   ) :: this
  integer,                     intent(in   ) :: unit, vlist(:)
  character(len=*),            intent(in   ) :: iotype
  integer,                     intent(  out) :: iostat
  character(len=*),            intent(inout) :: iomsg

  integer :: n, m


  m = size(this%psiN)
  n = size(this%turns)
  WRITETXT(metadata_fmt("PSIN", "i0"), m)
  WRITETXT(metadata_fmt("TURNS", "i0"), n)
  WRITETXT(metadata_fmt("NSAMPLES", "i0"), this%nsamples)
  WRITETXT(ewd_fmt(m, vlist, .true.), this%psin)
  WRITETXT(ewd_fmt(m, vlist, .true.), this%turns)
  WRITETXT(iwm_fmt(n+1, vlist, .true.), this%counts(:,:,-1))
  WRITETXT(iwm_fmt(n+1, vlist), this%counts(:,:,1))

  end subroutine write_formatted
  !-----------------------------------------------------------------------------

end module flare_connection_length
