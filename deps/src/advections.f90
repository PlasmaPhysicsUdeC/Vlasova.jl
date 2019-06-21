! Compile with flags: -shared -fPIC
! Ex: gfortran -Wall -shared -fimplicit-none -fopenmp -fPIC test.f90 -o test.so

! 1D advections

! subroutine velocity_advection1d( Nx, Nvx2p1, coef, Ex, ux, F)
!   implicit none
!   integer, intent(in):: Nx, Nvx2p1
!   real*8, intent(in)::  Ex(Nx), ux(Nvx2p1)
!   complex*16, intent(in):: coef
!   complex*16, intent(inout):: F(Nx, Nvx2p1)

!   integer:: i, j
!   real*8:: u

!   !$OMP PARALLEL DO PRIVATE(i,j,u)
!   do j = 1, Nvx2p1
!      u = ux(j)
!      do i = 1, Nx
!         F(i, j) = F(i, j) * exp( coef * Ex(i) * u )
!      end do
!   end do

!   return
! end subroutine velocity_advection1d

subroutine velocity_advection1d( Nx, Nvx2p1, coef, Ex, ux, F)
  implicit none
  integer, intent(in):: Nx, Nvx2p1
  real*8, intent(in)::  Ex(Nx), ux(Nvx2p1)
  complex*16, intent(in):: coef
  complex*16, intent(inout):: F(Nx, Nvx2p1)

  integer:: i, j
  complex*16:: tmp(Nx), tmp2(Nx)
  
  tmp = exp( coef * ( ux(2) - ux(1) ) * Ex )
  tmp2(:) = 1.0d0                  ! whole Array assignment

  !TODO: Is it convenient to parallelize inner loops? 
  do j = 1, Nvx2p1
     !$OMP PARALLEL DO PRIVATE(i)
     do i = 1, Nx
        F(i, j) = F(i, j) * tmp2(i)
        tmp2(i) = tmp2(i) * tmp(i)
     end do
  end do

  return
end subroutine velocity_advection1d


subroutine space_advection1d( Nx2p1, Nvx, spaceShift, F)
  implicit none
  integer, intent(in):: Nx2p1, Nvx
  complex*16, intent(in):: spaceShift(Nx2p1, Nvx)
  complex*16, intent(inout):: F(Nx2p1, Nvx)

  integer:: i, j

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j = 1, Nvx
     do i = 1, Nx2p1
        F(i, j) = F(i, j) * spaceShift(i, j)
     end do
  end do

  return
end subroutine space_advection1d

! 2D advections

subroutine velocity_advection2d( Nx, Ny, Nvx2p1, Nvy, coef, Ex, Ey, ux, uy, F)
  implicit none
  integer, intent(in):: Nx, Ny, Nvx2p1, Nvy
  real*8, intent(in)::  Ex(Nx, Ny), Ey(Nx, Ny), ux(Nvx2p1), uy(Nvy)
  complex*16, intent(in):: coef
  complex*16, intent(inout):: F(Nx, Ny, Nvx2p1, Nvy)

  integer:: i, j, k ,l
  real*8:: u1, u2

  !$OMP PARALLEL DO PRIVATE(i,j,k,l,u1,u2)
  do l = 1, Nvy
     u2 = uy(l)
     do k = 1, Nvx2p1
        u1 = ux(k)
        do j = 1, Ny
           do i = 1, Nx
              F(i, j, k, l) = F(i, j, k, l) * exp( coef * ( Ex(i, j) * u1 + Ey(i, j) * u2 ) )
           end do
        end do
     end do
  end do

  return
end subroutine velocity_advection2d

subroutine space_advection2d( Nx2p1, Ny, Nvx, Nvy, spaceShift, F)
  implicit none
  integer, intent(in):: Nx2p1, Ny, Nvx, Nvy
  complex*16, intent(in):: spaceShift(Nx2p1, Ny, Nvx, Nvy)
  complex*16, intent(inout):: F(Nx2p1, Ny, Nvx, Nvy)

  integer:: i, j, k ,l

  !$OMP PARALLEL DO PRIVATE(i,j,k,l)
  do l = 1, Nvy
     do k = 1, Nvx
        do j = 1, Ny
           do i = 1, Nx2p1
              F(i, j, k ,l) = F(i, j, k ,l) * spaceShift(i, j, k, l)
           end do
        end do
     end do
  end do

  return
end subroutine space_advection2d
