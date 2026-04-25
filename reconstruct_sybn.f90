subroutine reconstruct_vof_sybn(aL,aR,aa,dx,cell_type,nfirst,nlast,ngc,Nmat)

  use compute_method_choosing,only:limiter,Ord_lim
  implicit none

  integer :: nfirst, nlast, ngc, Nmat

  double precision aL(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision aR(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision aa(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision dx(nfirst-ngc:nlast+ngc)
  integer cell_type(nfirst-ngc:nlast+ngc)

  integer :: i,j,k,m,n

  double precision slop,slopL,slopR,dxL,dxR,tmp

  i=nfirst-ngc
  aL(:,i) = aa(:,i)
  aR(:,i) = aa(:,i)
  
  i=nlast+ngc
  aL(:,i) = aa(:,i)
  aR(:,i) = aa(:,i)
  
  do i=nfirst-ngc+1,nlast+ngc-1
     if(cell_type(i).eq.1) then

        aL(:,i) = aa(:,i)
        aR(:,i) = aa(:,i)

     else
        
        dxL= (dx(i-1)+dx(i))/2.d0
        dxR= (dx(i)+dx(i+1))/2.d0
        
        m=2
        slopL= (aa(m,i)-aa(m,i-1))/dxL
        slopR= (aa(m,i+1)-aa(m,i))/dxR
        if(limiter.eq.0) then
           slop = 0.d0
        else if(limiter.eq.1) then
           if(slopL*slopR.le.0.d0) then
              slop = 0.d0
           else
              slop = sign(dmin1(dabs(slopL),dabs(slopR)),slopL)
           end if
        else if(limiter.eq.2) then
           if(slopL*slopR.le.0.d0) then
              slop = 0.d0
           else
              slop = 2.d0*slopL*dabs(slopR)/(dabs(slopL)+dabs(slopR)+1.d-20)
           end if
        else
           print *, " + this case of limiter is not included : "
           print *, " + limiter =", limiter
           pause
           stop
        end if
  
        tmp = dx(i)
        aL(m,i)= aa(m,i) -slop/2.d0*tmp
        aR(m,i)= aa(m,i) +slop/2.d0*tmp

        aL(m,i)= dmin1(dmax1(aL(m,i),1.d-11),1.d0-1.d-11)
        aR(m,i)= dmin1(dmax1(aR(m,i),1.d-11),1.d0-1.d-11)

        aL(1,i)= 1.d0-aL(2,i)
        aR(1,i)= 1.d0-aR(2,i)

     end if
     
  end do
    
end subroutine reconstruct_vof_sybn


subroutine reconstruct_sybn(aL,aR,aa,vof,dx,cell_type,nfirst,nlast,ngc,Nmat)

  use compute_method_choosing,only:limiter,Ord_lim
  implicit none

  integer :: nfirst, nlast, ngc, Nmat

  double precision aL(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision aR(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision aa(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision vof(1:Nmat,nfirst-ngc:nlast+ngc)
  double precision dx(nfirst-ngc:nlast+ngc)
  integer cell_type(nfirst-ngc:nlast+ngc)

  integer :: i,j,k,m,n

  double precision slop,slopL,slopR,dxL,dxR,tmp,theta

  i=nfirst-ngc
  aL(:,i) = aa(:,i)
  aR(:,i) = aa(:,i)

  i=nlast+ngc
  aL(:,i) = aa(:,i)
  aR(:,i) = aa(:,i)

  do i=nfirst-ngc+1,nlast+ngc-1

     dxL= (dx(i-1)+dx(i))/2.d0
     dxR= (dx(i)+dx(i+1))/2.d0

! for material 1
     m=1
     slopL= (aa(m,i)-aa(m,i-1))/dxL
     slopR= (aa(m,i+1)-aa(m,i))/dxR
     if(limiter.eq.0) then
        slop = 0.d0
     else if(limiter.eq.1) then
        if(slopL*slopR.le.0.d0) then
           slop = 0.d0
        else
           slop = sign(dmin1(dabs(slopL),dabs(slopR)),slopL)
        end if
     else if(limiter.eq.2) then
        if(slopL*slopR.le.0.d0) then
           slop = 0.d0
        else
           slop = 2.d0*slopL*dabs(slopR)/(dabs(slopL)+dabs(slopR)+1.d-20)
        end if
     else
        print *, " + this case of limiter is not included : "
        print *, " + limiter =", limiter
        pause
        stop
     end if

     tmp = dx(i)
     aL(m,i)= aa(m,i) -slop/2.d0*tmp
     aR(m,i)= aa(m,i) +slop/2.d0*tmp

! for material 2
     m=2
     if(cell_type(i).eq.2.and.cell_type(i-1).eq.2.and. cell_type(i+1).eq.2) then
        slopL= (aa(m,i)-aa(m,i-1))/dxL
        slopR= (aa(m,i+1)-aa(m,i))/dxR
        if(limiter.eq.0) then
           slop = 0.d0
        else if(limiter.eq.1) then
           if(slopL*slopR.le.0.d0) then
              slop = 0.d0
           else
              slop = sign(dmin1(dabs(slopL),dabs(slopR)),slopL)
           end if
        else if(limiter.eq.2) then
           if(slopL*slopR.le.0.d0) then
              slop = 0.d0
           else
              slop = 2.d0*slopL*dabs(slopR)/(dabs(slopL)+dabs(slopR)+1.d-20)
           end if
        else
           print *, " + this case of limiter is not included : "
           print *, " + limiter =", limiter
           pause
           stop
        end if

        tmp = dmin1( dmin1(vof(2,i-1),vof(2,i+1))/dmax1(vof(2,i),1.d-99), 1.d0)   

        if(tmp.le.Ord_lim(1)) then
           theta = 0.d0
        else if(tmp.ge.Ord_lim(2)) then
           theta = 1.d0
        else
           theta = (tmp-Ord_lim(1))/dmax1(Ord_lim(2)-Ord_lim(1),1.d-99)
        end if

        slop= theta*slop 

     else

        slop =0.d0

     end if

     tmp = dx(i)
     aL(m,i)= aa(m,i) -slop/2.d0*tmp
     aR(m,i)= aa(m,i) +slop/2.d0*tmp
     
  end do

end subroutine reconstruct_sybn


