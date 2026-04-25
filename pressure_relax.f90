subroutine pressure_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
     xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

  !  use Local_variables 

  implicit none

  integer :: ifirst0, ilast0, gcw0, ifirst1, ilast1, gcw1, ifirst2, ilast2, gcw2, Nmat

  double precision &
       xp(ifirst0-gcw0-1:ilast0+gcw0), &
       dx(ifirst0-gcw0:ilast0+gcw0), &
       yp(ifirst1-gcw1-1:ilast1+gcw1), &
       dy(ifirst1-gcw1:ilast1+gcw1), &
       zp(ifirst2-gcw2-1:ilast2+gcw2), &
       dz(ifirst2-gcw2:ilast2+gcw2) 

  integer &
       ke(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2),&
       ep_index(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2)

  double precision &
       fc(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       dens(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       iener(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       vx(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       vy(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       vz(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       pr(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2), &
       snd(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2),&      
       tempr(Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2),&           
       stress(12,Nmat,ifirst0-gcw0:ilast0+gcw0,&
       ifirst1-gcw1:ilast1+gcw1,&
       ifirst2-gcw2:ilast2+gcw2)

  double precision dt

  integer :: i,j,k,i1,j1,k1,m,n 

  integer, dimension(2) :: ke_loc
  double precision, dimension(2):: vof_loc,dens_loc,iener_loc,pr_loc,snd_loc,tempr_loc

  k=ifirst2
  do j=ifirst1,ilast1
     do i=ifirst0,ilast0

        ke_loc(1:2) = ke(1:2,i,j,k)
        vof_loc(1:2) = fc(1:2,i,j,k)
        dens_loc(1:2) = dens(1:2,i,j,k)
        iener_loc(1:2) = iener(1:2,i,j,k)
        pr_loc(1:2) = pr(1:2,i,j,k)
        snd_loc(1:2) = snd(1:2,i,j,k)
        tempr_loc(1:2) = tempr(1:2,i,j,k)


        call stiff_pressure_relax_u2ph(vof_loc,dens_loc,iener_loc,pr_loc,snd_loc,tempr_loc,ke_loc)

        fc(1:2,i,j,k) = vof_loc(1:2)
        dens(1:2,i,j,k) = dens_loc(1:2)
        iener(1:2,i,j,k) = iener_loc(1:2)
        pr(1:2,i,j,k) = pr_loc(1:2)
        snd(1:2,i,j,k) = snd_loc(1:2)
        tempr(1:2,i,j,k) = tempr_loc(1:2)

     end do
  end do

end subroutine pressure_relax

