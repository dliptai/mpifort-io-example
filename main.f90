program Main
  use mpi
  use utils
  use io_write
  use io_read

  implicit none

  integer :: ierr, type_mpi_array
  integer :: i,j,k, ranki
  integer, parameter :: header1 = 1
  integer :: header1_tmp
  real(8), parameter :: header2 = 1.0
  real(8) :: header2_tmp
  character(len=*), parameter :: filename = 'output.dat'
  real(8), allocatable, dimension(:,:,:) :: x1, x2
  real(8), external :: xfunc
  logical :: legacy = .true.
  integer :: fh, error

  call init_mpi

  is_global = 1 ; ie_global = 8
  js_global = 1 ; je_global = 8
  ks_global = 1 ; ke_global = 8

  call decompose_domain
  call create_type_mpi_array(type_mpi_array)

  ! allocate and initialise array
  allocate(x1(is:ie,js:je,ks:ke))
  allocate(x2,mold=x1)

  do i = is,ie
    do j = js,je
      do k = ks,ke
        x1(i,j,k) = xfunc(i,j,k)
        x2(i,j,k) = -xfunc(i,j,k)
      end do
    end do
  end do

  ! create empty file on master task
  if (myrank == 0) then
    open(newunit=fh, file=filename, status='replace', form='unformatted', action='write', access='stream')
    close(unit=fh)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call open_file_write(filename, fh)
  call write_header(fh, header1, header2, legacy)
  call write_array_mpi(fh, x1, type_mpi_array, legacy)
  call write_array_mpi(fh, x2, type_mpi_array, legacy)
  call close_file(fh)
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  ! put junk in array
  x1 = -1.0
  x2 = -1.0

  ! -------------- Now read back the same file -----------------------------------------

  ! call is_file_legacy(filename, legacy) ! this doesn't work now, because the first thing in the header is an integer...

  call open_file_read(filename, fh)
  call read_header(fh, header1_tmp, header2_tmp, legacy)
  call read_array_mpi(fh, x1, type_mpi_array, legacy)
  call read_array_mpi(fh, x2, type_mpi_array, legacy)
  call close_file(fh)

  if (myrank ==0 ) write(*,'(A)') '>>> Checking header...'
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (header1 - header1_tmp /= 0) then
    write(*, '(A,I0,A,I0)') 'rank ', myrank, ' Header1 mismatch: ', header1, header1_tmp
  else
    write(*, '(A,I0,A,I0)') 'rank ', myrank, ' Header1 matches!'
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  if (abs(header2 - header2_tmp) > 1e-14) then
    write(*, '(A,I0,A,F0.1,A,F0.1)') 'rank ', myrank, ' Header2 mismatch: ', header2, header2_tmp
  else
    write(*, '(A,I0,A,F0.1,A,F0.1)') 'rank ', myrank, ' Header2 matches!'
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  error = 0

  if (myrank ==0 ) write(*,'(A)') '>>> Checking arrays...'
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  do ranki = 0, nprocs-1
    if (myrank == ranki) then
      do i = is,ie
        do j = js,je
          do k = ks,ke
            if (abs(x1(i,j,k) - xfunc(i,j,k)) > 1e-14) then
              print*, 'Rank ', myrank, ' x1 has a problem at ', i, j, k, x1(i,j,k), abs(x1(i,j,k) - xfunc(i,j,k))
              error = error + 1
            end if
            if (abs(x2(i,j,k) + xfunc(i,j,k)) > 1e-14) then
              print*, 'Rank ', myrank, ' x2 has a problem at ', i, j, k, x2(i,j,k), abs(x2(i,j,k) + xfunc(i,j,k))
              error = error + 1
            end if
          end do
        end do
      end do
    endif
    call mpi_barrier(MPI_COMM_WORLD, ierr)
  end do

  call mpi_barrier(MPI_COMM_WORLD, ierr)

  if (error == 0) then
    print*, 'Rank ', myrank, ' All tests passed'
  else
    print*, 'Rank ', myrank, ' N tests failed:', error
  end if

  call MPI_FINALIZE(ierr)

end program Main

real(8) function xfunc(i, j, k)
  implicit none
  integer, intent(in) :: i, j, k

  xfunc = real(100*i + 10*j + k, kind=8)

end function xfunc
