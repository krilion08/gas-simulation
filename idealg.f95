program ideal_gas
  use omp_lib
  implicit none
  integer, parameter :: natoms=2000
  real(kind=8), allocatable, dimension(:,:) :: r, v
  real(kind=8) :: start, finish, dt
  real(kind=8), parameter ::  eps=0.5
  real(kind=8), dimension(3) :: box_size=10.0d0
  integer :: axis, i, j, k, flag

  
  allocate (r(3,natoms), v(3,natoms))
  call random_number(r)
  call random_number(v)
  v=v-0.5
  r=(r+0.5)*5
  dt=0.15
  flag=0

  write(*,*) v(1,2), v(3,3), r(1,3), r(4,3)
  call write_xyz(r, flag, natoms, box_size, 0*dt)
  do i=1,50000 
    r=r+v*dt
    flag=1
    call write_xyz(r, flag, natoms, box_size, i*dt)
    if (any((r(:,:) .le. eps) .or. any(r(:,:) .ge. box_size(1)-eps))) then
      do j=1,natoms

        do k=1,3
          if ((r(k,j) .le. eps) .or. (r(k,j) .ge. box_size(k)-eps)) then
          call wall (v, j, k)
          endif
        enddo

      enddo
    endif 
  enddo
  deallocate(v,r)
end program



subroutine wall(v, particle, axis)
  implicit none
  real(kind=8), dimension(3,1) :: v
  integer :: axis, particle
  v(axis, particle)=-v(axis, particle)
end subroutine

subroutine write_xyz(r, flag, natoms, box_size, t)
  implicit none
  real(kind=8), dimension(3,natoms) :: r
  real(kind=8) :: t
  real(kind=8), dimension(3) :: box_size
  integer :: i, flag, natoms

  open(1, file = 'trajectory.xyz', action = 'write', position = 'append', status = 'unknown')
  write(1,101) natoms
  if (flag .eq. 0) then
    write(1,*) "simulation of ", natoms, " non interracting molecules in a box with dimensions of ", &
    box_size(1), "x", box_size(2),"x", box_size(3), ". t = ", t
  else
    write(1,*) "t = ", t
  endif
  do i=1,natoms
    write(1,103) "Xe", r(:,i)
  enddo
  close (1)

!xyz format1
101 format(I5)
!xyz format2
103 format(A2, 2X, F7.4, X, F7.4, X, F7.4)
end subroutine



