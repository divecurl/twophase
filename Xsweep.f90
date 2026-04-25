subroutine Xsweep(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
        xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)
    
  use compute_method_choosing,only: bound, Coord
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
  
  index_direction = 1

  is_boundary_1d = 1
  bound_1d(1) = bound(1)
  bound_1d(2) = bound(2)
  
  call boundary(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
        xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,1,0,0)
  
  k=ifirst2
  do j=ifirst1,ilast1
     y_panal = j
     z_panal = k

     xp_1d(-gcw0) = xp(ifirst0-gcw0-1)
     
     do i= ifirst0-gcw0,ilast0+gcw0
        i1= i-ifirst0+1
        xp_1d(i1) = xp(i)
        dx_1d(i1) = dx(i)
     end do     
     
     do i=ifirst0-gcw0,ilast0+gcw0
        i1= i-ifirst0+1

        ke_1d(1:Nmat,i1)= ke(1:Nmat,i,j,k)
        vof_1d(1:Nmat,i1)= fc(1:Nmat,i,j,k)

        dd_1d(1:Nmat,i1)= dens(1:Nmat,i,j,k)
        ee_1d(1:Nmat,i1)= iener(1:Nmat,i,j,k)
        pp_1d(1:Nmat,i1)= pr(1:Nmat,i,j,k)
        cc_1d(1:Nmat,i1)= snd(1:Nmat,i,j,k)
        tt_1d(1:Nmat,i1)= tempr(1:Nmat,i,j,k)
        vx_1d(1:Nmat,i1)= vx(1:Nmat,i,j,k)
        vy_1d(1:Nmat,i1)= vy(1:Nmat,i,j,k)
        vz_1d(1:Nmat,i1)= vz(1:Nmat,i,j,k)
        
     end do

     call one_dim_deltaU_u2ph(dt,1,ilast0-ifirst0+1,gcw0,index_direction)         

     do i=ifirst0-gcw0,ilast0+gcw0
        i1= i-ifirst0+1

        D_vof(1:Nmat,i,j,k) = D_vof(1:Nmat,i,j,k)+D_vof_1d(1:Nmat,i1)
        D_mas(1:Nmat,i,j,k) = D_mas(1:Nmat,i,j,k)+D_mas_1d(1:Nmat,i1)
        D_xmo(1:Nmat,i,j,k) = D_xmo(1:Nmat,i,j,k)+D_xmo_1d(1:Nmat,i1)
        D_ymo(1:Nmat,i,j,k) = D_ymo(1:Nmat,i,j,k)+D_ymo_1d(1:Nmat,i1)
        D_zmo(1:Nmat,i,j,k) = D_zmo(1:Nmat,i,j,k)+D_zmo_1d(1:Nmat,i1)
        D_eng(1:Nmat,i,j,k) = D_eng(1:Nmat,i,j,k)+D_eng_1d(1:Nmat,i1)
        
     end do

  end do
   
end subroutine Xsweep
