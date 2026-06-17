!===============================================================================
! Table: 2D array of real values
!
!
! Constructors:
!    table(rows, columns, lbounds=[1,1])
!                               allocate new table of given shape and set values to zero
!
!    make(values, lbounds=[1,1])
!                               construct table from 2D array
!
!    load(filename, columns, transposed=.false., lbounds=[1,1])
!                               load table from data file
!
!
! Components:
!    values			elements in table
!
!
! Type-bound procedures:
!    broadcast()		broadcast table to all processes
!    allreduce()		sum data values over all processes (it is the users responsibility to ensure that table sizes are compatible on all processes)
!    write(iu)			write table to unit
!    rows()			return number of rows in table
!    columns()			return number of columns in table
!    sort(column)		sort table with respect to values in selected column
!    binary_search(column, T)	return row index where T matches (sorted) column values
!    reverse_rows()
!    reverse_columns()
!    append_rows(n)		append n rows to table
!    append_columns(m)		append m columns to table
!    prepend_rows(n)		prepend n rows to table
!    prepend_columns(m)		prepend m columns to table
!    remove_row(i)
!    remove_column(j)
!    remove_rows(n)		remove n rows from end of table
!    remove_columns(m)		remove m columns from end of table
!    resize_rows(n)		resize table to n rows by appending or removing rows
!    resize_columns(m)		resize table to m columns by appending or removing columns
!    resize(n,m)		resize table to n rows and m columns
!===============================================================================
module moose_table
  use iso_fortran_env
  use moose_txtio
  implicit none
  private


  type, extends(txtio), public :: table
     real(real64), dimension(:,:), allocatable :: values

     contains
     ! broadcast table
     procedure :: broadcast

     ! finalize table
     procedure :: free

     ! sum data values over all processes
     procedure :: allreduce

     ! write table
     procedure :: write_selection
     procedure :: write_formatted

     ! return number of rows/columns in table
     procedure :: rows
     procedure :: columns

     ! advanced data manipulation
     procedure :: sort
     procedure :: binary_search
     procedure :: reverse_rows
     procedure :: reverse_columns

     ! manage table shape
     procedure :: append_rows
     procedure :: append_columns
     procedure :: prepend_rows
     procedure :: prepend_columns

     procedure :: remove_row
     procedure :: remove_column

     procedure :: remove_rows
     procedure :: remove_columns

     procedure :: resize_rows
     procedure :: resize_columns
     procedure :: resize
  end type table


  interface table
     procedure new
     procedure make
     procedure load
  end interface



  interface resize_array
     procedure :: resize_int1, resize_int2, resize_int3
     procedure :: resize_real1, resize_real2, resize_real3
  end interface resize_array


  public :: &
     resize_array

  contains
  !-----------------------------------------------------------------------------


! constructors:
  !-----------------------------------------------------------------------------
  function new(rows, columns, lbounds) result(T)
  !
  ! allocate new table of given shape and set values to zero
  !
  ! **parameters**:
  !
  ! :rows:        number of rows in table
  !
  ! :columns:     number of colums in table
  !
  ! **optional parameters:**
  !
  ! :lbounds:     lower bounds for values array (default: [1, 1])
  !
  integer, intent(in) :: rows, columns
  integer, intent(in), optional :: lbounds(2)
  type(table)         :: T

  integer :: lb(2), ub(2)


  call init_txtio(T, "table")
  if (rows <= 0  .or.  columns <= 0) then
     return
  endif
  lb = 1;  if (present(lbounds)) lb = lbounds
  ub = lb + [rows, columns] - 1

  allocate (T%values(lb(1):ub(1), lb(2):ub(2)), source=0.d0)

  end function new
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function make(values, lbounds) result(T)
  !
  ! construct table from 2D array
  !
  real(real64), intent(in) :: values(:,:)
  integer,      intent(in), optional :: lbounds(2)
  type(table)              :: T


  T = table(size(values,1), size(values,2), lbounds)
  T%values = values

  end function make
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  function load(filename, columns, transposed, lbounds) result(T)
  !
  ! load table from data file.
  !
  ! **parameters:**
  !
  ! :filename:     name of the data file
  !
  ! **optional parameters:**
  !
  ! :columns:      number of data columns to expect in file, exceeding columns are ignored (default: derived from first data line)
  !
  ! :transposed:   table is transposed after loading (default = .false.)
  !
  ! :lbounds:      lower bounds for values array (default: [1, 1])
  !
  character(len=*), intent(in) :: filename
  integer,          intent(in), optional :: columns
  logical,          intent(in), optional :: transposed
  integer,          intent(in), optional :: lbounds(2)
  type(table)                  :: T

  integer, parameter :: iu        = 99  ! unit number associated with data file
  integer, parameter :: ncol_max  = 256 ! max. number of columns supported
  integer, parameter :: col_width = 30  ! max. number of characters per column

  integer, parameter :: n = ncol_max * col_width

  character(len=n)   :: line
  real(real64)       :: r(ncol_max)
  logical            :: ex, lT
  integer            :: i, idim, iline, ios, nrow, ncol


  ! check if data file exists
  inquire(file=filename, exist=ex)
  if (.not.ex) then
     write (6, 9001) trim(filename);   stop
  endif
  open  (iu, file=filename)
 9001 format("error in table constructor load: data file ", a, " does not exist!")


  ! set default values for optional arguments
  ncol = 0;         if (present(columns))    ncol = columns
  lT   = .false.;   if (present(transposed)) lT   = transposed


  ! obtain number of columns from 1st data line (if not set by user)
  if (ncol == 0) then
     ! get 1st data line
     do
        read  (iu, 1000, iostat=ios) line
        if (ios /= 0) then
           write (6, 9002) filename;   stop
        endif
        if (line(1:1) /= '#') exit
     enddo

     ! process 1st data line to find number of columns
     do ncol=0,ncol_max-1
        read  (line, *, iostat=ios) r(1:ncol+1)
        if (ios /= 0) exit
     enddo
     rewind(iu)
  endif
  if (ncol <= 0) then
     write (6, 9003) ncol;   stop
  endif
 9002 format("error in table constructor load: unexpected end of data file ", a, "!")
 9003 format("error in table constructor load: invalid number of columns = ", i0)


  ! parse data file to get number of data lines -> rows in table
  nrow = 0
  parse_loop: do
     read  (iu, 1000, end=1001) line
     if (line(1:1) /= '#'  .and.  line /= "") nrow = nrow + 1
  enddo parse_loop
 1000 format(a)
 1001 rewind(iu)


  ! initialize table
  if (lT) then
     T = table(ncol, nrow, lbounds)
     idim = 2
  else
     T = table(nrow, ncol, lbounds)
     idim = 1
  endif


  ! load actual data
  iline = 0
  do i=lbound(T%values, idim),ubound(T%values, idim)
     ! get next data line
     do
        read  (iu, 1000) line
        iline = iline + 1
        if (line(1:1) /= '#' .and. line /= "") exit
     enddo

     ! process data line
     if (lT) then
        read (line, *, iostat=ios) T%values(:,i)
     else
        read (line, *, iostat=ios) T%values(i,:)
     endif
     if (ios /= 0) then
        write (6, 9004) ncol, iline, trim(filename)
        stop
     endif
  enddo
  close (iu)
 9004 format("error in table constructor load: expecting ", i0, " values in line ", i0, " of data file ", a, "!")

  end function load
  !-----------------------------------------------------------------------------


! type-bound procedures:
  !-----------------------------------------------------------------------------
  subroutine broadcast(this)
  use moose_mpi
  class(table), intent(inout) :: this


  call this%txtio_broadcast()
  call proc(0)%broadcast_allocatable(this%values)

  end subroutine broadcast
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine free(this)
  class(table), intent(inout) :: this


  deallocate (this%values)
  call this%txtio_free()

  end subroutine free
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine allreduce(this)
  !
  ! sum data values over all processes
  ! NOTE: it is the users responsibility to ensure that table sizes are compatible on all processes!
  !
  use moose_mpi
  class(table), intent(inout) :: this


  call moose_mpi_sum(this%values)

  end subroutine allreduce
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine write_selection(this, iu, vlist, iostat, iomsg, rows, columns, transposed)
  class(table),     intent(in   ) :: this
  integer,          intent(in   ) :: iu, vlist(:)
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg
  integer,          intent(in   ), optional :: rows(2), columns(2)
  logical,          intent(in   ), optional :: transposed

  integer :: i, i1, i2, j, j1, j2, n
  logical :: T


  ! set range for output rows and columns
  i1 = lbound(this%values, 1)
  i2 = ubound(this%values, 1)
  j1 = lbound(this%values, 2)
  j2 = ubound(this%values, 2)
  if (present(rows)) then
     i1 = max(i1, rows(1))
     i2 = min(i2, rows(2))
  endif
  if (present(columns)) then
     j1 = max(j1, columns(1))
     j2 = min(j2, columns(2))
  endif


  T = .false.;   if (present(transposed)) T = transposed
  if (T) then
     ! write transposed output
     n = i2-i1+1
     write (iu, ewd_fmt(n, vlist), iostat=iostat, iomsg=iomsg) (this%values(i1:i2,j), j=j1,j2)
  else
     ! write regular output
     n = j2-j1+1
     write (iu, ewd_fmt(n, vlist), iostat=iostat, iomsg=iomsg) (this%values(i,j1:j2), i=i1,i2)
  endif

  end subroutine write_selection
  !-----------------------------------------------------------------------------
  subroutine write_formatted(this, unit, iotype, vlist, iostat, iomsg)
  class(table),     intent(in   ) :: this
  integer,          intent(in   ) :: unit, vlist(:)
  character(len=*), intent(in   ) :: iotype
  integer,          intent(  out) :: iostat
  character(len=*), intent(inout) :: iomsg


  call this%write_selection(unit, vlist, iostat, iomsg)

  end subroutine write_formatted
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function rows(this)
  class(table), intent(in) :: this
  integer                  :: rows
  rows = size(this%values, 1) ! number of rows in table
  end function rows
  !-----------------------------------------------------------------------------
  pure function columns(this)
  class(table), intent(in) :: this
  integer                  :: columns
  columns = size(this%values, 2) ! number of columns in table
  end function columns
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine sort(this, jpivot, axis, irange)
  !
  ! sort table rows (axis = 1) or columns (axis = 2) with respect to values in
  ! jpivot-th column or row, respectively
  !
  use moose_algorithms, only: quicksort
  class(table), intent(inout) :: this
  integer,      intent(in   ) :: jpivot
  integer,      intent(in   ), optional :: axis, irange(2)

  integer :: i(2), axis_


  ! set axis
  axis_ = 1;   if (present(axis)) axis_ = axis

  ! set range
  i = [1, size(this%values, axis_)]
  if (present(irange)) i = 1 + irange - lbound(this%values, axis_)

  call quicksort(this%values, i(1), i(2), jpivot, axis_)

  end subroutine sort
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure function binary_search(this, column, T) result(m)
  use moose_algorithms, binary_search_array => binary_search_R
  class(table), intent(in) :: this
  integer,      intent(in) :: column
  real(real64), intent(in) :: T
  integer                  :: m


  m = binary_search_array(this%values(:,column), T) - 1 + lbound(this%values,1)

  end function binary_search
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine reverse_rows(this, first_row, last_row)
  class(table), intent(inout) :: this
  integer,      intent(in   ), optional :: first_row, last_row

  real(real64) :: tmp(size(this%values,2))
  integer :: i, i1, i2, j


  i1 = lbound(this%values,1);   if (present(first_row)) i1 = max(i1, first_row)
  i2 = ubound(this%values,1);   if (present(last_row))  i2 = min(i2, last_row)
  do i=i1,i1+(i2-i1+1)/2-1
     j   = i2-i+i1
     tmp = this%values(i,:)
     this%values(i,:) = this%values(j,:)
     this%values(j,:) = tmp
  enddo

  end subroutine reverse_rows
  !-----------------------------------------------------------------------------
  pure subroutine reverse_columns(this, first_column, last_column)
  class(table), intent(inout) :: this
  integer,      intent(in   ), optional :: first_column, last_column

  real(real64) :: tmp(size(this%values,1))
  integer :: i, i1, i2, j


  i1 = lbound(this%values,2);   if (present(first_column)) i1 = max(i1, first_column)
  i2 = ubound(this%values,2);   if (present(last_column))  i2 = min(i2, last_column)
  do i=i1,i1+(i2-i1+1)/2-1
     j   = i2-i+i1
     tmp = this%values(:,i)
     this%values(:,i) = this%values(:,j)
     this%values(:,j) = tmp
  enddo

  end subroutine reverse_columns
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine append_rows(this, rows)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: rows

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: l1, l2, u1, u2


  if (rows <= 0) return
  l1 = lbound(this%values,1)
  l2 = lbound(this%values,2)
  u1 = ubound(this%values,1)
  u2 = ubound(this%values,2)
  allocate (tmp(l1:u1+rows, l2:u2), source=0.d0)
  tmp(l1:u1, l2:u2) = this%values
  call move_alloc(tmp, this%values)

  end subroutine append_rows
  !-----------------------------------------------------------------------------
  pure subroutine append_columns(this, ncols)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: ncols

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: l1, l2, u1, u2


  if (ncols <= 0) return
  l1 = lbound(this%values,1)
  l2 = lbound(this%values,2)
  u1 = ubound(this%values,1)
  u2 = ubound(this%values,2)
  allocate (tmp(l1:u1, l2:u2+ncols), source=0.d0)
  tmp(l1:u1, l2:u2) = this%values
  call move_alloc(tmp, this%values)

  end subroutine append_columns
  !-----------------------------------------------------------------------------
  pure subroutine prepend_rows(this, rows)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: rows

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: l1, l2, u1, u2


  if (rows <= 0) return
  l1 = lbound(this%values,1)
  l2 = lbound(this%values,2)
  u1 = ubound(this%values,1)
  u2 = ubound(this%values,2)
  allocate (tmp(l1-rows:u1, l2:u2), source=0.d0)
  tmp(l1:u1, l2:u2) = this%values
  call move_alloc(tmp, this%values)

  end subroutine prepend_rows
  !-----------------------------------------------------------------------------
  pure subroutine prepend_columns(this, ncols)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: ncols

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: l1, l2, u1, u2


  if (ncols <= 0) return
  l1 = lbound(this%values,1)
  l2 = lbound(this%values,2)
  u1 = ubound(this%values,1)
  u2 = ubound(this%values,2)
  allocate (tmp(l1:u1, l2-ncols:u2), source=0.d0)
  tmp(l1:u1, l2:u2) = this%values
  call move_alloc(tmp, this%values)

  end subroutine prepend_columns
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine remove_row(this, i)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: i

  integer :: i1, i2


  i1 = lbound(this%values,1)
  i2 = ubound(this%values,1)
  if (i < i1  .or.  i > i2) return


  if (i < i2) then
     this%values(i:i2-1,:) = this%values(i+1:i2,:)
  endif
  call this%remove_rows(1)

  end subroutine remove_row
  !-----------------------------------------------------------------------------
  pure subroutine remove_column(this, j)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: j

  integer :: j1, j2


  j1 = lbound(this%values,2)
  j2 = ubound(this%values,2)
  if (j < j1  .or.  j > j2) return


  if (j < j2) then
     this%values(:,j:j2-1) = this%values(:,j+1:j2)
  endif
  call this%remove_columns(1)

  end subroutine remove_column
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine remove_elements(this, n, iside, idimension)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: n, iside, idimension

  real(real64), dimension(:,:), allocatable :: tmp
  integer :: l(2), u(2), i


  if (n <= 0) return
  if (iside /= 1  .and.  iside /= -1) return
  if (idimension < 1  .or.  idimension > 2) return
  do i=1,2
     if (i == idimension) then
        if (iside == 1) then
           l(i) = lbound(this%values, i)
           u(i) = max(ubound(this%values, i)-n, l(i)-1)
        else
           u(i) = ubound(this%values, i)
           l(i) = min(lbound(this%values, i)+n, u(i)+1)
        endif

     else
        l(i) = lbound(this%values, i)
        u(i) = ubound(this%values, i)
     endif
  enddo
  allocate (tmp(l(1):u(1), l(2):u(2)), source=0.d0)
  tmp = this%values(l(1):u(1), l(2):u(2))
  call move_alloc(tmp, this%values)

  end subroutine remove_elements
  !-----------------------------------------------------------------------------
  pure subroutine remove_rows(this, nrows, iside)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: nrows
  integer,      intent(in), optional :: iside

  integer :: i


  i = 1;   if (present(iside)) i = iside
  call remove_elements(this, nrows, i, 1)

  end subroutine remove_rows
  !-----------------------------------------------------------------------------
  pure subroutine remove_columns(this, ncols, iside)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: ncols
  integer,      intent(in), optional :: iside

  integer :: i


  i = 1;   if (present(iside)) i = iside
  call remove_elements(this, ncols, i, 2)

  end subroutine remove_columns
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  pure subroutine resize_rows(this, rows, lrow)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: rows
  integer,      intent(in), optional :: lrow

  integer :: i1, i2, n, lb, ub


  n  = size(this%values, 1)
  lb = lbound(this%values, 1)
  ub = ubound(this%values, 1)
  i1 = lb;   if (present(lrow)) i1 = lrow
  i2 = i1 + rows - 1

  ! adjust upper boundary of table
  if (i2 < ub) then
     call this%remove_rows(ub-i2)
  else
     call this%append_rows(i2-ub)
  endif

  ! adjust lower boundary of table (in case lrow is set)
  if (i1 > lb) then
     call this%remove_rows(i1-lb, iside=-1)
  else
     call this%prepend_rows(lb-i1)
  endif

  end subroutine resize_rows
  !-----------------------------------------------------------------------------
  pure subroutine resize_columns(this, columns, lcol)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: columns
  integer,      intent(in), optional :: lcol

  integer :: j1, j2, m, lb, ub


  m  = size(this%values, 2)
  lb = lbound(this%values, 2)
  ub = ubound(this%values, 2)
  j1 = lb;   if (present(lcol)) j1 = lcol
  j2 = j1 + columns - 1

  ! adjust upper boundary of table
  if (j2 < ub) then
     call this%remove_columns(ub-j2)
  else
     call this%append_columns(j2-ub)
  endif

  ! adjust lower boundary of table (in case lcol is set)
  if (j1 > lb) then
     call this%remove_columns(j1-lb, iside=-1)
  else
     call this%prepend_columns(lb-j1)
  endif

  end subroutine resize_columns
  !-----------------------------------------------------------------------------
  pure subroutine resize(this, rows, columns, lrow, lcol)
  class(table), intent(inout) :: this
  integer,      intent(in)    :: rows, columns
  integer,      intent(in), optional :: lrow, lcol


  call this%resize_rows(rows, lrow)
  call this%resize_columns(columns, lcol)

  end subroutine resize
  !-----------------------------------------------------------------------------


! module procedures:
  !-----------------------------------------------------------------------------
  subroutine resize_int1(x, n, lb, ub, source)
  integer, allocatable, intent(inout) :: x(:)
  integer,              intent(in   ), optional :: n, lb, ub, source

  integer, allocatable :: tmp(:)
  integer :: k0, k0move, k1, k1move


  ! new lower bound
  k0 = lbound(x, 1)
  if (present(lb)) k0 = lb

  ! new upper bound
  k1 = ubound(x, 1)
  if (present(n)) then
     k1 = k0 + n - 1
  elseif (present(ub)) then
     k1 = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0:k1), source = source)
  else
     allocate (tmp(k0:k1), source = 0)
  endif

  k0move = max(k0, lbound(x, 1))
  k1move = min(k1, ubound(x, 1))
  tmp(k0move:k1move) = x(k0move:k1move)
  call move_alloc(tmp, x)

  end subroutine resize_int1
  !-----------------------------------------------------------------------------
  subroutine resize_int2(x, dim, n, lb, ub, source)
  integer, allocatable, intent(inout) :: x(:,:)
  integer,              intent(in   ) :: dim
  integer,              intent(in   ), optional :: n, lb, ub, source

  integer, allocatable :: tmp(:,:)
  integer :: k0(2), k0move, k1(2), k1move


  ! new lower bound
  k0 = lbound(x)
  if (present(lb)) k0(dim) = lb

  ! new upper bound
  k1 = ubound(x)
  if (present(n)) then
     k1(dim) = k0(dim) + n - 1
  elseif (present(ub)) then
     k1(dim) = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0(1):k1(1), k0(2):k1(2)), source = source)
  else
     allocate (tmp(k0(1):k1(1), k0(2):k1(2)), source = 0)
  endif

  k0move = max(k0(dim), lbound(x, dim))
  k1move = min(k1(dim), ubound(x, dim))
  select case (dim)
  case (1)
     tmp(k0move:k1move,:) = x(k0move:k1move,:)
  case (2)
     tmp(:,k0move:k1move) = x(:,k0move:k1move)
  end select
  call move_alloc(tmp, x)

  end subroutine resize_int2
  !-----------------------------------------------------------------------------
  subroutine resize_int3(x, dim, n, lb, ub, source)
  integer, allocatable, intent(inout) :: x(:,:,:)
  integer,              intent(in   ) :: dim
  integer,              intent(in   ), optional :: n, lb, ub, source

  integer, allocatable :: tmp(:,:,:)
  integer :: k0(3), k0move, k1(3), k1move


  ! new lower bound
  k0 = lbound(x)
  if (present(lb)) k0(dim) = lb

  ! new upper bound
  k1 = ubound(x)
  if (present(n)) then
     k1(dim) = k0(dim) + n - 1
  elseif (present(ub)) then
     k1(dim) = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0(1):k1(1), k0(2):k1(2), k0(3):k1(3)), source = source)
  else
     allocate (tmp(k0(1):k1(1), k0(2):k1(2), k0(3):k1(3)), source = 0)
  endif

  k0move = max(k0(dim), lbound(x, dim))
  k1move = min(k1(dim), ubound(x, dim))
  select case (dim)
  case (1)
     tmp(k0move:k1move,:,:) = x(k0move:k1move,:,:)
  case (2)
     tmp(:,k0move:k1move,:) = x(:,k0move:k1move,:)
  case (3)
     tmp(:,:,k0move:k1move) = x(:,:,k0move:k1move)
  end select
  call move_alloc(tmp, x)

  end subroutine resize_int3
  !-----------------------------------------------------------------------------


  !-----------------------------------------------------------------------------
  subroutine resize_real1(x, n, lb, ub, source)
  real(real64), allocatable, intent(inout) :: x(:)
  integer,                   intent(in   ), optional :: n, lb, ub
  real(real64),              intent(in   ), optional :: source

  real(real64), allocatable :: tmp(:)
  integer :: k0, k0move, k1, k1move


  ! new lower bound
  k0 = lbound(x, 1)
  if (present(lb)) k0 = lb

  ! new upper bound
  k1 = ubound(x, 1)
  if (present(n)) then
     k1 = k0 + n - 1
  elseif (present(ub)) then
     k1 = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0:k1), source = source)
  else
     allocate (tmp(k0:k1), source = 0.d0)
  endif

  k0move = max(k0, lbound(x, 1))
  k1move = min(k1, ubound(x, 1))
  tmp(k0move:k1move) = x(k0move:k1move)
  call move_alloc(tmp, x)

  end subroutine resize_real1
  !-----------------------------------------------------------------------------
  subroutine resize_real2(x, dim, n, lb, ub, source)
  real(real64), allocatable, intent(inout) :: x(:,:)
  integer,                   intent(in   ) :: dim
  integer,                   intent(in   ), optional :: n, lb, ub
  real(real64),              intent(in   ), optional :: source

  real(real64), allocatable :: tmp(:,:)
  integer :: k0(2), k0move, k1(2), k1move


  ! new lower bound
  k0 = lbound(x)
  if (present(lb)) k0(dim) = lb

  ! new upper bound
  k1 = ubound(x)
  if (present(n)) then
     k1(dim) = k0(dim) + n - 1
  elseif (present(ub)) then
     k1(dim) = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0(1):k1(1), k0(2):k1(2)), source = source)
  else
     allocate (tmp(k0(1):k1(1), k0(2):k1(2)), source = 0.d0)
  endif

  k0move = max(k0(dim), lbound(x, dim))
  k1move = min(k1(dim), ubound(x, dim))
  select case (dim)
  case (1)
     tmp(k0move:k1move,:) = x(k0move:k1move,:)
  case (2)
     tmp(:,k0move:k1move) = x(:,k0move:k1move)
  end select
  call move_alloc(tmp, x)

  end subroutine resize_real2
  !-----------------------------------------------------------------------------
  subroutine resize_real3(x, dim, n, lb, ub, source)
  real(real64), allocatable, intent(inout) :: x(:,:,:)
  integer,                   intent(in   ) :: dim
  integer,                   intent(in   ), optional :: n, lb, ub
  real(real64),              intent(in   ), optional :: source

  real(real64), allocatable :: tmp(:,:,:)
  integer :: k0(3), k0move, k1(3), k1move


  ! new lower bound
  k0 = lbound(x)
  if (present(lb)) k0(dim) = lb

  ! new upper bound
  k1 = ubound(x)
  if (present(n)) then
     k1(dim) = k0(dim) + n - 1
  elseif (present(ub)) then
     k1(dim) = ub
  endif


  ! copy data and move allocation
  if (present(source)) then
     allocate (tmp(k0(1):k1(1), k0(2):k1(2), k0(3):k1(3)), source = source)
  else
     allocate (tmp(k0(1):k1(1), k0(2):k1(2), k0(3):k1(3)), source = 0.d0)
  endif

  k0move = max(k0(dim), lbound(x, dim))
  k1move = min(k1(dim), ubound(x, dim))
  select case (dim)
  case (1)
     tmp(k0move:k1move,:,:) = x(k0move:k1move,:,:)
  case (2)
     tmp(:,k0move:k1move,:) = x(:,k0move:k1move,:)
  case (3)
     tmp(:,:,k0move:k1move) = x(:,:,k0move:k1move)
  end select
  call move_alloc(tmp, x)

  end subroutine resize_real3
  !-----------------------------------------------------------------------------

end module moose_table
