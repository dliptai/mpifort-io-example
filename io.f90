module io
  use vars
  use utils
  implicit none

contains
  subroutine write_hormone_dump(filename)
    use io_write
    character(len=*) :: filename
    integer :: fh

    call open_file_write(filename, fh)
    call write_fake_recordmarker(fh)
    call write_header(fh, tn, time)
    call write_fake_recordmarker(fh)
    call write_fake_recordmarker(fh)
    call write_array_mpi(fh, d, type_mpi_array)
    call write_array_mpi(fh, v1, type_mpi_array)
    call write_array_mpi(fh, v2, type_mpi_array)
    call write_array_mpi(fh, v3, type_mpi_array)
    call write_array_mpi(fh, e, type_mpi_array)
    call write_fake_recordmarker(fh)
    call close_file(fh)

  end subroutine write_hormone_dump

  subroutine read_hormone_dump(filename)
    use io_read
    character(len=*) :: filename
    integer :: fh

    call open_file_read(filename, fh)
    offset = 0
    call read_fake_recordmarker(fh)
    call read_header(fh, tn, time)
    call read_fake_recordmarker(fh)
    call read_fake_recordmarker(fh)
    call read_array_mpi(fh, d, type_mpi_array)
    call read_array_mpi(fh, v1, type_mpi_array)
    call read_array_mpi(fh, v2, type_mpi_array)
    call read_array_mpi(fh, v3, type_mpi_array)
    call read_array_mpi(fh, e, type_mpi_array)
    call read_fake_recordmarker(fh)
    call close_file(fh)

  end subroutine read_hormone_dump
end module io
