module flare_control
  use iso_fortran_env
  use moose_mpi
  use moose_dict
  implicit none


  ! control paramters for screen output
  logical :: verbose = .false., very_verbose = .false.
  logical :: report = .false.


  ! control parameter for troubleshooting
  logical :: diagnostic_mode = .false.


  ! numerical parameters for separatrix2d
  real(real64) :: &
     separatrix2d_step_size = 1.d-2, &
     separatrix2d_offset    = 1.d-3, &
     separatrix2d_fX        = 0.95d0, &
     separatrix2d_epsabs    = 1.d-8, &
     separatrix2d_alpha     = 1.d0

  integer :: &
     separatrix2d_nmax      = 16


  ! statistics
  integer, parameter :: &
     BFIELD_EVAL        = 1, &
     JACOBIAN_EVAL      = 2, &
     FIELDLINE_EVOLVE   = 3, &
     BOUNDARY_INTERSECT = 4
  integer(kind=selected_int_kind(16)) :: counter(4) = 0


  ! path prefix for model database (directory tree)
  type(dict) :: database

  contains
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine flare_init(standalone, greeting)
  use moose_mpi
  use moose_configparser
  use moose_utils, only: expanduser
  use flare_version
  logical, intent(in), optional :: standalone, greeting

  logical, save :: initialized = .false.

  type(configparser) :: cp
  type(cp_section), pointer :: database_section
  type(dict_item), pointer :: item
  logical :: standalone_, greeting_


  if (initialized) return
  initialized = .true.

  ! initialize MPI
  standalone_ = .true.;   if (present(standalone)) standalone_ = standalone
  if (standalone_) then
     call moose_mpi_init()
  else
     call moose_mpi_setup()
  endif

  ! initialize screen output parameters
  greeting_ = .true.;   if (present(greeting)) greeting_ = greeting
  report = rank == 0
  if (report .and. greeting_) then
     print *
     print 1000
     print *, "FLARE - The Field Line Analysis and Reconstruction Environment (" // &
        achar(27) // "[94m" // version // achar(27) // "[0m)"

     if (nproc > 1) then
        print *
        print 1001, nproc
     endif
  endif
 1000 format(80("="))
 1001 format(3x,"- Parallelization with ",i0," processes")


  ! read user configuration
  if (rank == 0) then
     ! set default value
     call database%set("default", expanduser("~/DATABASE/flare"))

     ! parse configuration file
     cp = configparser()
     call cp%read(expanduser("~/.flare"))
     database_section => cp%find_section("database")
     if (associated(database_section)) then
        item => database_section%options%first_item()
        do
           if (.not.associated(item)) exit
           call database%set(item%key, expanduser(item%val))
           item => item%next
        enddo
     endif
  endif

  end subroutine flare_init
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine flare_finalize()
  use moose_mpi


  ! if (verbose) print statistics ...

  call moose_mpi_finalize()

  end subroutine flare_finalize
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make_idir(direction, reference) result(idir)
  !
  ! return sign (-1/1) for direction given by string
  !
  character(len=*), intent(in) :: direction
  integer,          intent(in) :: reference
  integer                      :: idir

  character(len=32) :: direction_


  ! verify valid reference
  if (abs(reference) /= 1) then
     print 9000, reference
     stop
  endif
 9000 format("ERROR: invalid reference '",i0,"'")


  ! expand abbreviations
  direction_ = direction
  select case(direction)
  case("fwd")
     direction_ = "forward"
  case("bwd")
     direction_ = "backward"
  case("cw")
     direction_ = "clockwise"
  case("ccw")
     direction_ = "counter-clockwise"
  end select


  ! set idir
  select case(direction_)
  case("forward")
     idir =  1
  case("backward")
     idir = -1
  case("clockwise")
     idir = -reference
  case("counter-clockwise")
     idir =  reference
  case default
     print 9001, trim(direction_)
     stop
  end select
  if (report) print 1000, trim(direction_)
 1000 format(3x,"- Direction: ",a)
 9001 format("ERROR: invalid direction '",a,"'")

  end function make_idir
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine increment_counter(key)
  use moose_error
  integer, intent(in) :: key


  if (key <= 0  .or.  key > size(counter)) then
     call INDEX_ERROR("key", [1, size(counter)], "increment_counter")
  endif
  counter(key) = counter(key) + 1

  end subroutine increment_counter
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine progress_bar(iteration, maximum)
  use moose_utils, moose_progress_bar => progress_bar
  integer, intent(in) :: iteration, maximum


  if (report .and. .not.verbose) call moose_progress_bar(iteration, maximum)

  end subroutine progress_bar
  !-----------------------------------------------------------------------------
  subroutine finalize_progress_bar()
  use moose_utils


  call wait_for_all_procs(report .and. .not.verbose)

  end subroutine finalize_progress_bar
  !-----------------------------------------------------------------------------

end module flare_control
