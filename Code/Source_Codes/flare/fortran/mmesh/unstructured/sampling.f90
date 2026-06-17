module flare_mmesh_unstructured_sampling
  use iso_fortran_env
  use flare_mmesh_unstructured_mmesh
  implicit none
  private


  type, public :: source
     ! total source strength and number of samples
     real(real64) :: total_src, w0, wreduce
     integer :: nsamples

     ! iside = 0:    volume source
     !         1-4:  surface source
     integer, allocatable :: itube(:)
     integer :: iside, nsrc

     ! internal reference to mesh
     type(mmesh), pointer :: mesh

     contains
     procedure :: sample, boundary_event, free
  end type source


  public :: &
     surface_source

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function surface_source(mesh, itag, iside, total_src, nsamples) result(this)
  class(mmesh), target, intent(in) :: mesh
  integer,              intent(in) :: itag, iside, nsamples
  real(real64),         intent(in) :: total_src
  type(source)                     :: this

  integer, parameter :: &
     COUNT_STAGE = 1, &
     DEFINE_STAGE = 2

  integer :: itube, istage


  this%total_src = total_src
  this%nsamples = nsamples
  this%mesh => mesh
  this%w0 = this%total_src / this%nsamples
  this%wreduce = 0.1d0 * this%w0

  this%iside = iside
  do istage=COUNT_STAGE,DEFINE_STAGE
     this%nsrc = 0
     do itube=0,mesh%ntubes-1
        if (mesh%next_tube(1, iside, itube) == 0  .and. &
            mod(mesh%next_tube(2, iside, itube), 10) == itag) then
           if (istage == DEFINE_STAGE) this%itube(this%nsrc) = itube
           this%nsrc = this%nsrc + 1
        endif
     enddo
     if (istage == COUNT_STAGE) allocate (this%itube(0:this%nsrc-1))
  enddo

  end function surface_source
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine sample(this, c, w)
  class(source), intent(inout) :: this
  type(mcoords), intent(  out) :: c
  real(real64),  intent(  out) :: w

  real(real64) :: rsrc, rphi
  integer :: iphi_tube(0:1), isrc


  call random_number(rsrc)
  isrc = int(rsrc * this%nsrc)
  c%itube = this%itube(isrc)
  w = this%w0

  call random_number(rphi)
  iphi_tube = this%mesh%iphi_zone(:, this%mesh%izone_tube(c%itube))
  c%iphi = int(iphi_tube(0) + rphi * (iphi_tube(1) - iphi_tube(0)))
  call random_number(c%t)

  select case(this%iside)
  case (0)
     call random_number(c%xi)

  case (1)
     call random_number(c%xi(1))
     c%xi(2) = -0.999999d0
  end select

  end subroutine sample
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine boundary_event(this, w, pload, wload)
  !
  ! handle boundary event with loss probability *p*
  !
  ! input:
  !    w          old/new particle weight
  !    pload      loss probability at boundary
  !
  ! output:
  !    wload      contribution to boundary load
  !
  class(source), intent(in   ) :: this
  real(real64),  intent(inout) :: w
  real(real64),  intent(in   ) :: pload
  real(real64),  intent(  out) :: wload

  real(real64) :: r


  ! continuous reduction of weight
  if (w > this%wreduce) then
     wload = w * pload

  ! binary event: continue with same weight or stop
  else
     call random_number(r)
     if (r < pload) then
        wload = w
        w = 0.d0
        return
     else
        wload = 0.d0
     endif
  endif
  w = w - wload

  end subroutine boundary_event
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(source), intent(inout) :: this


  deallocate (this%itube)

  end subroutine free
  !-----------------------------------------------------------------------------

end module flare_mmesh_unstructured_sampling
