subroutine mesher(ifirst0,ilast0,gcw0,&
  ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2,&
  xp,dx,yp,dy,zp,dz,xyz_lo,xyz_hi)

  implicit none

  integer :: ifirst0, ilast0, gcw0, ifirst1, ilast1, gcw1, ifirst2, ilast2, gcw2

  double precision &
       xp(ifirst0-gcw0-1:ilast0+gcw0), &
       dx(ifirst0-gcw0:ilast0+gcw0), &
       yp(ifirst1-gcw1-1:ilast1+gcw1), &
       dy(ifirst1-gcw1:ilast1+gcw1), &
       zp(ifirst2-gcw2-1:ilast2+gcw2), &
       dz(ifirst2-gcw2:ilast2+gcw2) 

  double  precision :: xyz_lo(0:2), xyz_hi(0:2)
  integer :: i,j,k
  double precision :: tmp

  print *
  print *, "+ generaing mesh ..."
  
  tmp = (xyz_hi(0) - xyz_lo(0))/real(ilast0-ifirst0+1)
  do i = ifirst0-gcw0-1, ilast0+gcw0
     xp(i)= xyz_lo(0) + real(i- (ifirst0 - 1))*tmp
  end do
  dx = tmp

  tmp = (xyz_hi(1) - xyz_lo(1))/real(ilast1-ifirst1+1)
  do i = ifirst1-gcw1-1, ilast1+gcw1
     yp(i)= xyz_lo(1) + real(i- (ifirst1 - 1))*tmp
  end do
  dy = tmp

  tmp = (xyz_hi(2) - xyz_lo(2))/real(ilast2-ifirst2+1)
  do i = ifirst2-gcw2-1, ilast2+gcw2
     zp(i)= xyz_lo(2) + real(i- (ifirst2 - 1))*tmp
  end do
  dz = tmp
       
end subroutine mesher
