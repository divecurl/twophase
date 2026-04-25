subroutine velocity_relax(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
        xp,dx,yp,dy,zp,dz,ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,dt)
    
  use compute_method_choosing,only: Job_index, bound, C_d, Coord, Diam_par
  
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

  integer :: i,j,k  

  double precision tmp, tmp1, tmp2

  double precision vof_gas, Cd, Cd_plus, D_par

  double precision dvel_abs, dvel_abs0, dvel_norm(2), dvel(2), mas(2), vel(2,2), vel_bck(2,2), tmom(2)

  double precision Re, vof_cr, vof_sss, mu_g
  
  integer index_case
  
  index_case = 2
  

  Cd= C_d
  D_par = Diam_par

  mu_g= 18.d-10

  vof_cr= 0.63d0
  vof_sss= 1.d0/(1.d0+0.02d0*(1.d0-vof_cr)/vof_cr)
  
  k=ifirst2
  do j=ifirst1,ilast1
     do i=ifirst0, ilast0
         
        if(ke(2,i,j,k).ge.1) then

           vel_bck(1,1) = vx(1,i,j,k)
           vel_bck(2,1) = vy(1,i,j,k)
           vel_bck(1,2) = vx(2,i,j,k)
           vel_bck(2,2) = vy(2,i,j,k)

           mas(1) = fc(1,i,j,k)*dens(1,i,j,k)
           mas(2) = fc(2,i,j,k)*dens(2,i,j,k)

           tmom(:)= mas(1)*vel_bck(:,1)+mas(2)*vel_bck(:,2)

           dvel(:)= vel_bck(:,1)-vel_bck(:,2)
           dvel_abs0= dsqrt(dvel(1)**2.d0+dvel(2)**2.d0)
           dvel_norm(:)=dvel(:)/dmax1(dvel_abs0, 1.d-99) 

           vof_gas= fc(1,i,j,k)

           if(index_case.eq.1) then

              Re= fc(1,i,j,k)*dens(1,i,j,k)*dvel_abs0*D_par/mu_g
              Re=dmax1(Re,1.d-9)
              if(fc(2,i,j,k).gt.vof_cr) then
                 Cd= 150.d0*fc(2,i,j,k)/Re+1.75d0
              else if(fc(2,i,j,k).gt.1.d0-vof_sss) then
                 Cd= 150.d0*fc(2,i,j,k)/Re+1.75d0*((1.d0-vof_cr)*fc(2,i,j,k)/vof_cr/fc(1,i,j,k))**0.45d0
              else
                 Cd= 150d0*fc(2,i,j,k)/Re+0.3d0
              end if


              Cd_plus = Cd*dens(1,i,j,k)*fc(2,i,j,k)/D_par  

              
           else if(index_case.eq.2) then


              tmp1= 1.75d0
              tmp2= 0.75d0*0.44d0*fc(1,i,j,k)**(-1.65d0)
              
              tmp = fc(2,i,j,k)

              tmp = atan(150d0*1.75d0*(0.2d0-tmp))/3.1415926d0 + 0.5d0
              
              Cd= (1.d0-tmp)*tmp1 + tmp*tmp2  
              
              Cd_plus = Cd*dens(1,i,j,k)*fc(2,i,j,k)/D_par   

               
           else
               
              Cd= C_d ! from input file
              Cd_plus = Cd*dens(1,i,j,k)*fc(2,i,j,k)/D_par

           end if
           
           Cd_plus= Cd_plus*(1.d0/mas(1)+1.d0/mas(2))

           dvel_abs= 1.d0/dmax1(1.d0/dvel_abs0+Cd_plus*dt, 1.d-99)

           dvel(:)= dvel_abs*dvel_norm(:)

           vel(:,1) = (dvel(:)*mas(2)+tmom(:))/(mas(1)+mas(2))

           vel(:,2) = vel(:,1)- dvel(:)

           vx(1,i,j,k)= vel(1,1)
           vy(1,i,j,k)= vel(2,1)
           vx(2,i,j,k)= vel(1,2)
           vy(2,i,j,k)= vel(2,2)

           iener(1,i,j,k) = iener(1,i,j,k) + (mas(2)/(mas(1)+mas(2)))**2.d0/2.d0*(dvel_abs0**2.d0- dvel_abs**2.d0)
           iener(2,i,j,k) = iener(2,i,j,k) + (mas(1)/(mas(1)+mas(2)))**2.d0/2.d0*(dvel_abs0**2.d0- dvel_abs**2.d0)

           call eos_press(dens(1,i,j,k),iener(1,i,j,k),pr(1,i,j,k),snd(1,i,j,k),tempr(1,i,j,k),tmp,ke(1,i,j,k))
           call eos_press(dens(2,i,j,k),iener(2,i,j,k),pr(2,i,j,k),snd(2,i,j,k),tempr(2,i,j,k),tmp,ke(2,i,j,k))

        end if

     end do
  end do
  
end subroutine velocity_relax
