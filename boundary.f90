subroutine boundary(ifirst0,ilast0,gcw0,ifirst1,ilast1,gcw1,ifirst2,ilast2,gcw2, Nmat,&
     xp,dx,yp,dy,zp,dz, ke,fc,dens,iener,vx,vy,vz,pr,snd,tempr,ep_index,stress,id_x,id_y,id_z)

  use compute_method_choosing,only: bound

  implicit none

  integer :: ifirst0, ilast0, gcw0, ifirst1, ilast1, gcw1, ifirst2, ilast2, gcw2, Nmat

  integer :: id_x,id_y,id_z

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

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc x - direction 

  if(id_x.eq.1) then 

     ! x direction: left      
     if(bound(1).eq.1) then  
        k=ifirst2
        do j=ifirst1,ilast1
           do i=ifirst0-gcw0,ifirst0-1
              !i1=2*ifirst0-1-i
              i1=ifirst0
              
              ke(:,i,j,k)= ke(:,i1,j,k)
              ep_index(:,i,j,k)= ep_index(:,i1,j,k)
              fc(:,i,j,k)= fc(:,i1,j,k)
              dens(:,i,j,k)= dens(:,i1,j,k)
              iener(:,i,j,k)= iener(:,i1,j,k)
              vx(:,i,j,k)= vx(:,i1,j,k)
              vy(:,i,j,k)= vy(:,i1,j,k)
              vz(:,i,j,k)= vz(:,i1,j,k)
              pr(:,i,j,k)= pr(:,i1,j,k)
              snd(:,i,j,k)= snd(:,i1,j,k)
              tempr(:,i,j,k)= tempr(:,i1,j,k)
              stress(:,:,i,j,k)= stress(:,:,i1,j,k)
!              stress(2,:,i,j,k)= stress(2,:,i1,j,k)         
           end do
        end do

     else if(bound(1).eq.-1) then    

        k=ifirst2
        do j=ifirst1,ilast1
           do i=ifirst0-gcw0,ifirst0-1
              i1=2*ifirst0-1-i
              ke(:,i,j,k)= ke(:,i1,j,k)
              ep_index(:,i,j,k)= ep_index(:,i1,j,k)
              fc(:,i,j,k)= fc(:,i1,j,k)
              dens(:,i,j,k)= dens(:,i1,j,k)
              iener(:,i,j,k)= iener(:,i1,j,k)
              vx(:,i,j,k)= - vx(:,i1,j,k)
              vy(:,i,j,k)= vy(:,i1,j,k)
              vz(:,i,j,k)= vz(:,i1,j,k)
              pr(:,i,j,k)= pr(:,i1,j,k)
              snd(:,i,j,k)= snd(:,i1,j,k)
              tempr(:,i,j,k)= tempr(:,i1,j,k)
              stress(:,:,i,j,k)= stress(:,:,i1,j,k)
              stress(2,:,i,j,k)= - stress(2,:,i1,j,k)         
           end do
        end do

     endif

     !  x-direction: right 

     if(bound(2).eq.1) then  !   

        k=ifirst2
        do j=ifirst1,ilast1
           do i=ilast0+1,ilast0+gcw0
!              i1=2*ilast0+1-i
              i1=ilast0
              
              ke(:,i,j,k)= ke(:,i1,j,k)
              ep_index(:,i,j,k)= ep_index(:,i1,j,k)
              fc(:,i,j,k)= fc(:,i1,j,k)
              dens(:,i,j,k)= dens(:,i1,j,k)
              iener(:,i,j,k)= iener(:,i1,j,k)
              vx(:,i,j,k)= vx(:,i1,j,k)
              vy(:,i,j,k)= vy(:,i1,j,k)
              vz(:,i,j,k)= vz(:,i1,j,k)
              pr(:,i,j,k)= pr(:,i1,j,k)
              snd(:,i,j,k)= snd(:,i1,j,k)
              tempr(:,i,j,k)= tempr(:,i1,j,k)
              stress(:,:,i,j,k)= stress(:,:,i1,j,k)
           end do
        end do

     else if(bound(2).eq.-1) then    

        k=ifirst2
        do j=ifirst1,ilast1
           do i=ilast0+1,ilast0+gcw0
              i1=2*ilast0+1-i
              ke(:,i,j,k)= ke(:,i1,j,k)
              ep_index(:,i,j,k)= ep_index(:,i1,j,k)
              fc(:,i,j,k)= fc(:,i1,j,k)
              dens(:,i,j,k)= dens(:,i1,j,k)
              iener(:,i,j,k)= iener(:,i1,j,k)
              vx(:,i,j,k)= - vx(:,i1,j,k)
              vy(:,i,j,k)= vy(:,i1,j,k)
              vz(:,i,j,k)= vz(:,i1,j,k)
              pr(:,i,j,k)= pr(:,i1,j,k)
              snd(:,i,j,k)= snd(:,i1,j,k)
              tempr(:,i,j,k)= tempr(:,i1,j,k)
              stress(:,:,i,j,k)= stress(:,:,i1,j,k)
              stress(2,:,i,j,k)= - stress(2,:,i1,j,k)         
           end do
        end do

     endif

  end if



  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc y - direction 

  if(id_y.eq.1) then 

     ! y direction : blow      
     if(bound(3).eq.1) then  
        k=ifirst2
        do i=ifirst0,ilast0
           do j=ifirst1-gcw1,ifirst1-1

              j1=2*ifirst1-1-j
              ke(:,i,j,k)= ke(:,i,j1,k)
              ep_index(:,i,j,k)= ep_index(:,i,j1,k)
              fc(:,i,j,k)= fc(:,i,j1,k)
              dens(:,i,j,k)= dens(:,i,j1,k)
              iener(:,i,j,k)= iener(:,i,j1,k)
              vx(:,i,j,k)= vx(:,i,j1,k)
              vy(:,i,j,k)= vy(:,i,j1,k)
              vz(:,i,j,k)= vz(:,i,j1,k)
              pr(:,i,j,k)= pr(:,i,j1,k)
              snd(:,i,j,k)= snd(:,i,j1,k)
              tempr(:,i,j,k)= tempr(:,i,j1,k)
              stress(:,:,i,j,k)= stress(:,:,i,j1,k)

           end do
        end do

     else if(bound(3).eq.-1) then    

        k=ifirst2
        do i=ifirst0,ilast0
           do j=ifirst1-gcw1,ifirst1-1

              j1=2*ifirst1-1-j
              ke(:,i,j,k)= ke(:,i,j1,k)
              ep_index(:,i,j,k)= ep_index(:,i,j1,k)
              fc(:,i,j,k)= fc(:,i,j1,k)
              dens(:,i,j,k)= dens(:,i,j1,k)
              iener(:,i,j,k)= iener(:,i,j1,k)
              vx(:,i,j,k)= vx(:,i,j1,k)
              vy(:,i,j,k)= - vy(:,i,j1,k)
              vz(:,i,j,k)= vz(:,i,j1,k)
              pr(:,i,j,k)= pr(:,i,j1,k)
              snd(:,i,j,k)= snd(:,i,j1,k)
              tempr(:,i,j,k)= tempr(:,i,j1,k)
              stress(:,:,i,j,k)= stress(:,:,i,j1,k)
              stress(2,:,i,j,k)= - stress(2,:,i,j1,k)         

           end do
        end do

     endif

     !  y-direction: up 

     if(bound(4).eq.1) then  !   

        k=ifirst2
        do i=ifirst0,ilast0
           do j= ilast1+1,ilast1+gcw1

              j1=2*ilast1+1-j
              ke(:,i,j,k)= ke(:,i,j1,k)
              ep_index(:,i,j,k)= ep_index(:,i,j1,k)
              fc(:,i,j,k)= fc(:,i,j1,k)
              dens(:,i,j,k)= dens(:,i,j1,k)
              iener(:,i,j,k)= iener(:,i,j1,k)
              vx(:,i,j,k)= vx(:,i,j1,k)
              vy(:,i,j,k)= vy(:,i,j1,k)
              vz(:,i,j,k)= vz(:,i,j1,k)
              pr(:,i,j,k)= pr(:,i,j1,k)
              snd(:,i,j,k)= snd(:,i,j1,k)
              tempr(:,i,j,k)= tempr(:,i,j1,k)
              stress(:,:,i,j,k)= stress(:,:,i,j1,k)

           end do
        end do

     else if(bound(4).eq.-1) then    

        k=ifirst2
        do i=ifirst0,ilast0
           do j= ilast1+1,ilast1+gcw1

              j1=2*ilast1+1-j
              ke(:,i,j,k)= ke(:,i,j1,k)
              ep_index(:,i,j,k)= ep_index(:,i,j1,k)
              fc(:,i,j,k)= fc(:,i,j1,k)
              dens(:,i,j,k)= dens(:,i,j1,k)
              iener(:,i,j,k)= iener(:,i,j1,k)
              vx(:,i,j,k)= vx(:,i,j1,k)
              vy(:,i,j,k)= - vy(:,i,j1,k)
              vz(:,i,j,k)= vz(:,i,j1,k)
              pr(:,i,j,k)= pr(:,i,j1,k)
              snd(:,i,j,k)= snd(:,i,j1,k)
              tempr(:,i,j,k)= tempr(:,i,j1,k)
              stress(:,:,i,j,k)= stress(:,:,i,j1,k)
              stress(2,:,i,j,k)= - stress(2,:,i,j1,k)         
           end do
        end do

     endif

  end if


  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc z - direction 

  !need to be completed




end subroutine boundary
