program Main
  use mpi
  use vars
  use utils
  use io

  implicit none

  integer :: ierr
  integer :: i,j,k
  character(len=*), parameter :: filename = 'output.dat'
  real(8), external :: xfunc
  integer :: fh

  call init_mpi

  ! Set bounds of full domain
  is_global = 1 ; ie_global = 8
  js_global = 1 ; je_global = 8
  ks_global = 1 ; ke_global = 8

  call decompose_domain
  call create_type_mpi_array

  ! allocate and initialise arrays
  allocate(d(is:ie,js:je,ks:ke))
  allocate(v1,v2,v3,e,mold=d)
  do i = is,ie
    do j = js,je
      do k = ks,ke
        d(i,j,k)  = 1*xfunc(i,j,k)
        v1(i,j,k) = 2*xfunc(i,j,k)
        v2(i,j,k) = 3*xfunc(i,j,k)
        v3(i,j,k) = 4*xfunc(i,j,k)
        e(i,j,k)  = 5*xfunc(i,j,k)
      end do
    end do
  end do

  tn = 10
  time = 111.0

  ! create empty file on master task
  if (myrank == 0) then
    open(newunit=fh, file=filename, status='replace', form='unformatted', action='write', access='stream')
    close(unit=fh)
  end if
  call mpi_barrier(MPI_COMM_WORLD, ierr)

  call write_hormone_dump(filename)
  call read_hormone_dump(filename)

  if (tn /= 10) then
    print *, 'Error: tn /= 10'
  end if

  if ( abs(time - 111.0) > 1.e-14) then
    print *, 'Error: time /= 111.0'
  end if

  do i = is,ie
    do j = js,je
      do k = ks,ke
        if ( abs(d(i,j,k) - 1*xfunc(i,j,k)) > 1.e-14) then
          print *, 'Error: d /= 1*xfunc', 'at:', i,j,k
        end if
        if ( abs(v1(i,j,k) - 2*xfunc(i,j,k)) > 1.e-14) then
          print *, 'Error: v1 /= 2*xfunc', 'at:', i,j,k
        end if
        if ( abs(v2(i,j,k) - 3*xfunc(i,j,k)) > 1.e-14) then
          print *, 'Error: v2 /= 3*xfunc', 'at:', i,j,k
        end if
        if ( abs(v3(i,j,k) - 4*xfunc(i,j,k)) > 1.e-14) then
          print *, 'Error: v3 /= 4*xfunc', 'at:', i,j,k
        end if
        if ( abs(e(i,j,k) - 5*xfunc(i,j,k)) > 1.e-14) then
          print *, 'Error: e /= 5*xfunc', 'at:', i,j,k
        end if
      end do
    end do
  end do

  call MPI_FINALIZE(ierr)

end program Main

real(8) function xfunc(i, j, k)
  implicit none
  integer, intent(in) :: i, j, k

  xfunc = real(100*i + 10*j + k, kind=8)

end function xfunc
