module sli1

  use number
  use numerical

  contains

subroutine ferdi(l,m,wan,ngrid,stepvec,Rnl)
!
!JDY: 2021.2.22
!call Ferdis relation.
!
! Rnl(rvec) = int(dOmega wnR(rvec) Y^*lm(theta,phi)
!
! where Omega is solid angle.
!

  implicit none

!  complex(8),allocatable,intent(in) :: wan(:,:,:)
  real(8),allocatable,intent(in) :: wan(:,:,:)
  real(8),intent(in) :: stepvec(3,3)
  integer(4),intent(in) :: l,m,ngrid(3)

  complex(8),allocatable,intent(out) :: Rnl(:,:,:)

  integer(4) :: i,j,k,nx,ny,nz
  real(8) :: x,y,z,r,theta,phi
  complex(8) :: valZlm
  
  allocate(Rnl(ngrid(1),ngrid(2),ngrid(3)))

  nx = int(ngrid(1)/2)
  ny = int(ngrid(2)/2)
  nz = int(ngrid(3)/2)


  do i = -nx,nx
    do j = -ny,ny
      do k = -nz,nz
        x = dble(i)*stepvec(1,1) + dble(i)*stepvec(2,1) + dble(i)*stepvec(3,1)
        y = dble(j)*stepvec(1,2) + dble(j)*stepvec(2,2) + dble(j)*stepvec(3,2)
        z = dble(k)*stepvec(1,3) + dble(k)*stepvec(2,3) + dble(k)*stepvec(3,3)
        call cart2polar(x,y,z,r,theta,phi) 
        call Zlm(l,m,theta,phi,valZlm)
        Rnl(i+nx,j+ny,k+nz) = wan(i+nx,j+ny,k+nz)*valZlm
      end do
    end do
  end do

end subroutine ferdi

subroutine grid_vec(ngrid,stepvec,vecout)

  implicit none

  real(8),intent(in) :: stepvec(3,3)
  real(8) :: x,y,z
  real(8),allocatable :: tmpvec(:,:,:,:)
  real(8),allocatable,intent(out) :: vecout(:,:,:,:)
  integer(4),intent(in) :: ngrid(3)
  integer(4) :: i,j,k,nx,ny,nz,m1,m2,m3

!add one points which is stripped for PBC.
  do i = 1, 3
    ngrid(i) = ngrid(i) + 1
  end do 

!for even number grid, we must set -1 because points are 2n+1.
  if (mod(ngrid(1),2)==0) then
    m1 = -1
  else
    m1 = 0
  end if
  if (mod(ngrid(2),2)==0) then
    m2 = -1
  else
    m2 = 0
  end if
  if (mod(ngrid(3),2)==0) then
    m3 = -1
  else
    m3 = 0
  end if


  nx = int(ngrid(1)/2)
  ny = int(ngrid(2)/2)
  nz = int(ngrid(3)/2)

  allocate(vecout(ngrid(1),ngrid(2),ngrid(3),3))

  do i = -nx,nx + m1
    do j = -ny,ny + m2
      do k = -nz,nz + m3
        x = dble(i)*stepvec(1,1) + dble(i)*stepvec(2,1) + dble(i)*stepvec(3,1)
        y = dble(j)*stepvec(1,2) + dble(j)*stepvec(2,2) + dble(j)*stepvec(3,2)
        z = dble(k)*stepvec(1,3) + dble(k)*stepvec(2,3) + dble(k)*stepvec(3,3)
        vecout(i+nx+1,j+ny+1,k+nz+1,1) = x
        vecout(i+nx+1,j+ny+1,k+nz+1,2) = y
        vecout(i+nx*1,j+ny+1,k+nz+1,3) = z
      end do
    end do
  end do

end subroutine tran_vec

end module sli1
