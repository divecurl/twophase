subroutine cfl_1d_dt(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2,Nmat,xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,dt)
  
  use compute_method_choosing, only: cfl
  
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
       ifirst2-gcw2:ilast2+gcw2) 
     
  double precision dt
  double precision tmp, tmp1, tmp2, tmp3  
  double precision hh, v_max, c_max 
  integer :: i,j,k,n,m
  
  dt= 1.d8
  k=ifirst2
  j= ifirst1
  do i= ifirst0, ilast0
        
  v_max = 1.d-99
  c_max = 1.d-99
  do m= 1, Nmat
     if(ke(m,i,j,k).ge.1) then
        v_max = max(v_max, dabs(vx(m,i,j,k)) )
        c_max = max(c_max, snd(m,i,j,k))
     endif
  enddo
  
  tmp = cfl* dx(i)/(c_max+v_max)
  
  dt = min(dt,tmp)
  
  end do       
       
end subroutine cfl_1d_dt
    
    
subroutine cfl_2d_dt(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2,Nmat,xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,dt)
  
  use compute_method_choosing, only: cfl
  
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
       ifirst2-gcw2:ilast2+gcw2) 
     
  double precision dt
  double precision tmp, tmp1, tmp2, tmp3  
  double precision hh, vx_max, vy_max, c_max 
  integer :: i,j,k,n,m
  
  dt= 1.d8
  k=ifirst2
  do j= ifirst1, ilast1
  do i= ifirst0, ilast0
        
  vx_max = 1.d-99
  vy_max = 1.d-99  
  c_max = 1.d-99
  do m= 1, Nmat
     if(ke(m,i,j,k).ge.1) then
        vx_max = max(vx_max, dabs(vx(m,i,j,k)))
        vy_max = max(vy_max, dabs(vy(m,i,j,k)))
        c_max = max(c_max, snd(m,i,j,k))        
     endif
  end do
  
  tmp1 = dx(i)/(c_max+vx_max)
  tmp2 = dy(j)/(c_max+vy_max)
  
  dt = min(dt,tmp1, tmp2)
  
  end do
  end do       
       
  dt= cfl*dt
  
end subroutine cfl_2d_dt
