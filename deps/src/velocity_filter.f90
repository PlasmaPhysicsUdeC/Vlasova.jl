! Compile with flags: -shared -fPIC
! Ex: gfortran -Wall -shared -fimplicit-none -fopenmp -fPIC test.f90 -o test.so

subroutine velocity_filter1d(Nx, Nvx2p1, filter, F)
  implicit none
  integer, intent(in):: Nx, Nvx2p1
  real*8, intent(in)::  filter(Nvx2p1)
  complex*16, intent(inout):: F(Nx, Nvx2p1)

  integer:: i, j
  real*8:: vfilter

  !$OMP PARALLEL DO PRIVATE(I,J,vfilter)
  do j = 1, Nvx2p1
    vfilter = filter(j)
    do i = 1, Nx
      F(i, j) = F(i, j) * vfilter
    end do
  end do

  return
end subroutine velocity_filter1d

subroutine velocity_filter2d(Nx, Ny, Nvx2p1, Nvy, filter1, filter2, F)
  implicit none
  integer, intent(in):: Nx, Ny, Nvx2p1, Nvy
  real*8, intent(in)::  filter1(Nvx2p1), filter2(Nvy)
  complex*16, intent(inout):: F(Nx, Ny, Nvx2p1, Nvy)

  integer:: i, j, k, l
  real*8:: vfilter
  
  !$OMP PARALLEL DO PRIVATE(I,J,K,L,vfilter)
  do l = 1, Nvy
     do k = 1, Nvx2p1
        vfilter = filter1(k)*filter2(l)
        do j = 1, Ny
           do i = 1, Nx
              F(i, j, k, l) = F(i, j, k, l) * vfilter
           end do
        end do
     end do
  end do
  
  return
end subroutine velocity_filter2d
