module utils
  use mpi

  implicit none

  integer(kind=MPI_OFFSET_KIND) :: bytes_read = 0
  integer(kind=MPI_OFFSET_KIND) :: bytes_written = 0

  ! Assuming fortran record markers are the size of an integer
  integer, parameter :: record_marker = sizeof(int(0,kind=4))

  integer :: myrank, nprocs
  integer :: nx, ny, nz
  integer :: is, ie, js, je, ks, ke
  integer :: is_global, ie_global, js_global, je_global, ks_global, ke_global

  contains

  subroutine init_mpi()
    integer :: ierr

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

  end subroutine init_mpi

  subroutine decompose_domain()
    integer :: dims(3), mycoords(3), cart_comm, ierr
    logical :: periods(3)

    nx = ie_global - is_global + 1
    ny = je_global - js_global + 1
    nz = ke_global - ks_global + 1

    dims = [nprocs, 1, 1] ! Trivial decomposition, for now... TODO
    periods = [.true., .true., .true.] ! Always set to periodic and allow boundary conditions to override

    call MPI_CART_CREATE(MPI_COMM_WORLD, 3, dims, periods, .false., cart_comm, ierr)
    call MPI_CART_COORDS(cart_comm, myrank, 3, mycoords, ierr)

    ! Compute local domain, including offset if starting index is not 1
    is = mycoords(1) * (nx / dims(1)) + 1 - (is_global - 1)
    ie = (mycoords(1) + 1) * (nx / dims(1)) - (is_global - 1)
    js = mycoords(2) * (ny / dims(2)) + 1 - (js_global - 1)
    je = (mycoords(2) + 1) * (ny / dims(2)) - (js_global - 1)
    ks = mycoords(3) * (nz / dims(3)) + 1 - (ks_global - 1)
    ke = (mycoords(3) + 1) * (nz / dims(3)) - (ks_global - 1)

    ! Special treatment of final proc, in case nx/ny/nz is not divisible by nprocs
    if (mycoords(1) == dims(1) - 1) ie = nx - (is_global - 1)
    if (mycoords(2) == dims(2) - 1) je = ny - (js_global - 1)
    if (mycoords(3) == dims(3) - 1) ke = nz - (ks_global - 1)

    if (myrank == 0) then
      print*, 'Global domain: ', is_global, ie_global, js_global, je_global, ks_global, ke_global
    endif

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    print*, 'Rank ', myrank, ' has domain ', is, ie, js, je, ks, ke

  end subroutine decompose_domain

  subroutine set_array_view(array_view)
    integer, intent(out) :: array_view
    integer :: ierr, sizes(3), subsizes(3), starts(3)

    sizes = [nx,ny,nz]
    subsizes = [ie-is+1,je-js+1,ke-ks+1]
    starts = [is-1,js-1,ks-1]

    call mpi_type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, array_view, ierr)
    call mpi_type_commit(array_view, ierr)

  end subroutine set_array_view

  ! subroutine to try reading the header in sequential format. If it succeeds it should be a legacy file
  subroutine is_file_legacy(filename, legacy)
    character(len=*), intent(in) :: filename
    logical, intent(out) :: legacy
    real(8) :: header
    integer :: ierr, un

    ! Try reading the header in sequential format
    open(newunit=un, file=filename, status='old', form='unformatted', action='read', access='sequential')
    read(un, iostat=ierr) header
    close(un)

    if (ierr /= 0) then
      legacy = .false.
    else
      legacy = .true.
    end if

    if (myrank == 0) print *, 'Legacy file =', legacy

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine is_file_legacy

  subroutine write_header(filename, header, legacy)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: header(:)
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: un, ierr
    if (present(legacy)) islegacy = legacy

    if (myrank == 0) then

      if (islegacy) then
        open(newunit=un, file=filename, status='replace', form='unformatted', action='write', access='sequential')
      else
        open(newunit=un, file=filename, status='replace', form='unformatted', action='write', access='stream')
      end if

      write(un) header
      close(un)
    end if

    if (islegacy) then
      bytes_written = bytes_written + (2*record_marker + sizeof(header))
    else
      bytes_written = bytes_written + sizeof(header)
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine write_header

  ! a function to determine the total number of elements of a 3D array
  function total_elements(array) result(n)
    integer :: n
    real(8), intent(in) :: array(:,:,:)

    n = size(array,1) * size(array,2) * size(array,3)

  end function total_elements

  subroutine write_array_mpi(filename, buffer, array_view, legacy)
    character(len=*), intent(in) :: filename
    real(8), intent(in) :: buffer(:,:,:)
    integer, intent(in) :: array_view
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: file, status(MPI_STATUS_SIZE), ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: extent
    integer(kind=MPI_OFFSET_KIND) :: disp = 0
    if (present(legacy)) islegacy = legacy

    if (islegacy) call write_fake_recordmarker(filename)

    call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY + MPI_MODE_CREATE + MPI_MODE_APPEND, MPI_INFO_NULL, file, ierr)
    call mpi_file_get_position(file, disp, ierr)
    call mpi_file_set_view(file, disp, MPI_DOUBLE_PRECISION, array_view, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_write_all(file, buffer, total_elements(buffer), MPI_DOUBLE_PRECISION, status, ierr)
    call mpi_file_get_type_extent(file, array_view, extent, ierr)
    bytes_written = bytes_written + extent
    call mpi_file_close(file, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (islegacy) call write_fake_recordmarker(filename)

  end subroutine write_array_mpi

  subroutine write_fake_recordmarker(filename)
    character(len=*), intent(in) :: filename
    integer :: un, ierr
    integer(4) :: rm = -10

    if (myrank == 0) then
      open(newunit=un, file=filename, status='old', form='unformatted', action='write', access='stream', position='append')
      write(un) rm
      close(un)
    end if

    bytes_written = bytes_written + record_marker
    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine write_fake_recordmarker

  subroutine read_header(filename, header, legacy)
    character(len=*), intent(in) :: filename
    real(8), intent(out) :: header(:)
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: ierr, un
    if (present(legacy)) islegacy = legacy

    ! Put some junk in the header to make sure it's overwritten
    header = -1.0

    if (islegacy) then
      open(newunit=un, file=filename, status='old', form='unformatted', action='read', access='sequential')
    else
      open(newunit=un, file=filename, status='old', form='unformatted', action='read', access='stream')
    end if

    read(un) header
    close(un)

    if (islegacy) then
      bytes_read = bytes_read + (2*record_marker + sizeof(header))
    else
      bytes_read = bytes_read + sizeof(header)
    end if

    call mpi_barrier(MPI_COMM_WORLD, ierr)

  end subroutine read_header

  subroutine read_array_mpi(filename, buffer, array_view, legacy)
    character(len=*), intent(in) :: filename
    real(8), intent(out) :: buffer(:,:,:)
    integer, intent(in) :: array_view
    logical, intent(in), optional :: legacy
    logical :: islegacy = .false.
    integer :: file, status(MPI_STATUS_SIZE), ierr
    INTEGER(KIND=MPI_ADDRESS_KIND) :: extent
    if (present(legacy)) islegacy = legacy

    if (islegacy) bytes_read = bytes_read + record_marker

    call mpi_file_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, file, ierr)
    call mpi_file_set_view(file, bytes_read, MPI_DOUBLE_PRECISION, array_view, 'native', MPI_INFO_NULL, ierr)
    call mpi_file_read_all(file, buffer, total_elements(buffer), MPI_DOUBLE_PRECISION, status, ierr)
    call mpi_file_get_type_extent(file, array_view, extent, ierr)
    bytes_read = bytes_read + extent
    call mpi_file_close(file, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)

    if (islegacy) bytes_read = bytes_read + record_marker

  end subroutine read_array_mpi

end module utils
