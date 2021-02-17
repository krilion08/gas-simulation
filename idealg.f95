program ideal_gas
  use omp_lib
  implicit none
  real(kind=8) :: r(3,1)
  real(kind=8) :: v(3,1)=0
  real(kind=8) :: start, finish, dt
  real(kind=8), parameter :: box_size=10.d0, eps=0.5
  real(kind=8), dimension(3) :: wall_collision
  integer :: axis, i, j

 

  call random_number(r)
  call random_number(v)
  v=v-0.5
  r=abs(r+0.5)*5
  dt=0.25

  write(*,*) "R(10) = ", r(:,1), v(:,1)
  call write_xyz(r)
  do i=1,5000 
    r=r+v*dt
!    write(*,*) "R(10) = ", r(:,1), v(:,1)
    call write_xyz(r)
    if (any((r(:,:) .le. eps) .or. any(r(:,:) .ge. 10-eps))) then
      do j=1,3
        if ((r(j,1) .le. eps) .or. (r(j,1) .ge. 10-eps)) then
          call wall (v, 1, j)
          exit
        endif 
      enddo
    endif
  enddo
end program



subroutine wall(v, particle, axis)
  implicit none
  real(kind=8), dimension(3,1) :: v
  integer :: axis, particle
  v(axis, particle)=-v(axis, particle)
end subroutine

subroutine write_xyz(r)
  implicit none
  real(kind=8), dimension(3,1) :: r

  open(1, file = 'trajectory.xyz', action = 'write', position = 'append', status = 'unknown')
  write(1,101) size(r(1,:))
  write(1,*)
  write(1,103) "Xe", r(:,1)
  close (1)
!xyz format1
101 format(I4)
!xyz format2
103 format(A2, 4X, F7.4, X, F7.4, X, F7.4)
end subroutine



