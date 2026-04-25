subroutine NonSplit_FV2D(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
     xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

  use compute_method_choosing,only:  Dimen, Pre_Relax, Vel_Relax, Vof_lim_val
  use Local_variables 

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

  double precision tmp

  k=ifirst2
  do j=ifirst1,ilast1
     do i=ifirst0,ilast0

        do n=1,2
           vof_rk(n,i,j,k)= fc(n,i,j,k)
           mas_rk(n,i,j,k)= vof_rk(n,i,j,k)*dens(n,i,j,k)
           xmo_rk(n,i,j,k)= mas_rk(n,i,j,k)*vx(n,i,j,k)
           ymo_rk(n,i,j,k)= mas_rk(n,i,j,k)*vy(n,i,j,k)
           zmo_rk(n,i,j,k)= mas_rk(n,i,j,k)*vz(n,i,j,k)
           tmp = vx(n,i,j,k)**2.d0+vy(n,i,j,k)**2.d0+vz(n,i,j,k)**2.d0
           eng_rk(n,i,j,k)= mas_rk(n,i,j,k)*(iener(n,i,j,k)+tmp/2.d0)
        end do

        vof(1:2,i,j,k)= vof_rk(1:2,i,j,k)
        mas(1:2,i,j,k)= mas_rk(1:2,i,j,k)
        xmo(1:2,i,j,k)= xmo_rk(1:2,i,j,k)
        ymo(1:2,i,j,k)= ymo_rk(1:2,i,j,k)
        zmo(1:2,i,j,k)= zmo_rk(1:2,i,j,k)
        eng(1:2,i,j,k)= eng_rk(1:2,i,j,k)

     end do
  end do

  do n=1,2

     k=ifirst2
     do j=ifirst1,ilast1
        do i=ifirst0,ilast0
           D_vof(:,i,j,k)= 0.d0
           D_mas(:,i,j,k)= 0.d0
           D_xmo(:,i,j,k)= 0.d0
           D_ymo(:,i,j,k)= 0.d0
           D_zmo(:,i,j,k)= 0.d0
           D_eng(:,i,j,k)= 0.d0
        end do
     end do

     call Xsweep(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
          xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

     if(Dimen.eq.2) then
        call Ysweep(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
             xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)
     endif

     k=ifirst2
     do j=ifirst1,ilast1
        do i=ifirst0,ilast0

           vof(1:2,i,j,k) = vof(1:2,i,j,k) + D_vof(1:2,i,j,k)

           ! vof check
           tmp = vof(1,i,j,k) + vof(2,i,j,k)
           if(dabs(tmp-1.d0).ge.1.d-8) then
              print *, " + summation of volume fraction is not 1 :"
              print *, " + sum = ", tmp
              pause
              stop
           end if
           vof(1:2,i,j,k) = vof(1:2,i,j,k)/tmp
           if(vof(2,i,j,k).gt.Vof_lim_val(2) ) then
              print *, "+ inreasonable vof: "
              print *, "+  i,j,k =", i,j,k
              print *, "+ vof(1:2) = ", vof
              pause
              stop
           end if


           mas(1:2,i,j,k) = mas(1:2,i,j,k) + D_mas(1:2,i,j,k)
           xmo(1:2,i,j,k) = xmo(1:2,i,j,k) + D_xmo(1:2,i,j,k)
           ymo(1:2,i,j,k) = ymo(1:2,i,j,k) + D_ymo(1:2,i,j,k)
           zmo(1:2,i,j,k) = zmo(1:2,i,j,k) + D_zmo(1:2,i,j,k)
           eng(1:2,i,j,k) = eng(1:2,i,j,k) + D_eng(1:2,i,j,k)
        end do
     end do

     if(n.eq.1) then

        k=ifirst2
        do j=ifirst1,ilast1
           do i=ifirst0,ilast0
              fc(1:2,i,j,k) = vof(1:2,i,j,k)
              dens(1:2,i,j,k) = mas(1:2,i,j,k)/fc(1:2,i,j,k)
              vx(1:2,i,j,k) = xmo(1:2,i,j,k)/mas(1:2,i,j,k)

              vy(1:2,i,j,k) = ymo(1:2,i,j,k)/mas(1:2,i,j,k)

              vz(1:2,i,j,k) = zmo(1:2,i,j,k)/mas(1:2,i,j,k)

              do m=1,2
                 tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
                 iener(m,i,j,k) = eng(m,i,j,k)/mas(m,i,j,k)- tmp
                 call eos_press(dens(m,i,j,k),iener(m,i,j,k),pr(m,i,j,k),snd(m,i,j,k),tempr(m,i,j,k),tmp,ke(m,i,j,k))
              end do
           end do
        end do

        if(Vel_Relax.eq.2) then

           call velocity_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
                xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

           k=ifirst2
           do j=ifirst1,ilast1
              do i=ifirst0,ilast0
                 mas(1:2,i,j,k) = dens(1:2,i,j,k)*fc(1:2,i,j,k)
                 xmo(1:2,i,j,k) = mas(1:2,i,j,k)*vx(1:2,i,j,k)
                 ymo(1:2,i,j,k) = mas(1:2,i,j,k)*vy(1:2,i,j,k)
                 zmo(1:2,i,j,k) = mas(1:2,i,j,k)*vz(1:2,i,j,k)
                 do m=1,2
                    tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
                    eng(m,i,j,k) = mas(m,i,j,k)*(iener(m,i,j,k)+tmp)
                 end do
              end do
           end do

        end if

        if(Pre_Relax.eq.2) then

           call pressure_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
                xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

           k=ifirst2
           do j=ifirst1,ilast1
              do i=ifirst0,ilast0
                 mas(1:2,i,j,k) = dens(1:2,i,j,k)*fc(1:2,i,j,k)
                 xmo(1:2,i,j,k) = mas(1:2,i,j,k)*vx(1:2,i,j,k)
                 ymo(1:2,i,j,k) = mas(1:2,i,j,k)*vy(1:2,i,j,k)
                 zmo(1:2,i,j,k) = mas(1:2,i,j,k)*vz(1:2,i,j,k)
                 do m=1,2
                    tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
                    eng(m,i,j,k) = mas(m,i,j,k)*(iener(m,i,j,k)+tmp)
                 end do
              end do
           end do

        end if
     end if
  end do

!!! update 
  k=ifirst2
  do j=ifirst1,ilast1
     do i=ifirst0,ilast0
        vof(1:2,i,j,k) = (vof(1:2,i,j,k) + vof_rk(1:2,i,j,k))/2.d0

        ! vof check
        tmp = vof(1,i,j,k) + vof(2,i,j,k)
        if(dabs(tmp-1.d0).ge.1.d-8) then
           print *, " + summation of volume fraction is not 1 :"
           print *, " + sum = ", tmp
           pause
           stop
        end if
        vof(1:2,i,j,k) = vof(1:2,i,j,k)/tmp
        if(vof(2,i,j,k).gt.Vof_lim_val(2) ) then
           print *, "+ 2, inreasonable vof: "
           print *, "+  i,j,k =", i,j,k
           print *, "+ vof(1:2) = ", vof
           pause
           stop
        end if

        mas(1:2,i,j,k) = (mas(1:2,i,j,k) + mas_rk(1:2,i,j,k))/2.d0
        xmo(1:2,i,j,k) = (xmo(1:2,i,j,k) + xmo_rk(1:2,i,j,k))/2.d0
        ymo(1:2,i,j,k) = (ymo(1:2,i,j,k) + ymo_rk(1:2,i,j,k))/2.d0
        zmo(1:2,i,j,k) = (zmo(1:2,i,j,k) + zmo_rk(1:2,i,j,k))/2.d0
        eng(1:2,i,j,k) = (eng(1:2,i,j,k) + eng_rk(1:2,i,j,k))/2.d0

        fc(1:2,i,j,k) = vof(1:2,i,j,k)
        dens(1:2,i,j,k) = mas(1:2,i,j,k)/fc(1:2,i,j,k)
        vx(1:2,i,j,k) = xmo(1:2,i,j,k)/mas(1:2,i,j,k)
        vy(1:2,i,j,k) = ymo(1:2,i,j,k)/mas(1:2,i,j,k)
        vz(1:2,i,j,k) = zmo(1:2,i,j,k)/mas(1:2,i,j,k)
        do m=1,2
           tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
           iener(m,i,j,k) = eng(m,i,j,k)/mas(m,i,j,k)- tmp

           call eos_press(dens(m,i,j,k),iener(m,i,j,k),pr(m,i,j,k),snd(m,i,j,k),tempr(m,i,j,k),tmp,ke(m,i,j,k))
        end do

     end do
  end do


  if(Vel_Relax.ge.1) then
     call velocity_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
          xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

     k=ifirst2
     do j=ifirst1,ilast1
        do i=ifirst0,ilast0
           mas(1:2,i,j,k) = dens(1:2,i,j,k)*fc(1:2,i,j,k)
           xmo(1:2,i,j,k) = mas(1:2,i,j,k)*vx(1:2,i,j,k)
           ymo(1:2,i,j,k) = mas(1:2,i,j,k)*vy(1:2,i,j,k)
           zmo(1:2,i,j,k) = mas(1:2,i,j,k)*vz(1:2,i,j,k)
           do m=1,2
              tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
              eng(m,i,j,k) = mas(m,i,j,k)*(iener(m,i,j,k)+tmp)
           end do
        end do
     end do
  end if


  if(Pre_Relax.ge.1) then
     call pressure_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
          xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)

     k=ifirst2
     do j=ifirst1,ilast1
        do i=ifirst0,ilast0
           mas(1:2,i,j,k) = dens(1:2,i,j,k)*fc(1:2,i,j,k)
           xmo(1:2,i,j,k) = mas(1:2,i,j,k)*vx(1:2,i,j,k)
           ymo(1:2,i,j,k) = mas(1:2,i,j,k)*vy(1:2,i,j,k)
           zmo(1:2,i,j,k) = mas(1:2,i,j,k)*vz(1:2,i,j,k)
           do m=1,2
              tmp = (vx(m,i,j,k)**2.d0+vy(m,i,j,k)**2.d0+vz(m,i,j,k)**2.d0)/2.d0
              eng(m,i,j,k) = mas(m,i,j,k)*(iener(m,i,j,k)+tmp)
           end do
        end do
     end do
  end if

end subroutine NonSplit_FV2D_sybn
