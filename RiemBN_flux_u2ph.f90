subroutine RiemBN_flux_u2ph(flux_la,flux_eu,theL,dL,eL,pL,cL,vxL,vyL,vzL,keL,theR,dR,eR,pR,cR,vxR,vyR,vzR,keR)

  use compute_method_choosing,only: flux_solver,eos_stat
  implicit none

  double precision flux_la(0:5), flux_eu(0:6)

  double precision theL,dL,eL,pL,cL,vxL,vyL,vzL,theR,dR,eR,pR,cR,vxR,vyR,vzR


  integer :: keL, keR

  double precision av_L(2), av_R(2)

  double precision GamL, GamR

  double precision tmp, tmp1, tmp2
  integer :: itmp, n_step, eos_type, index_ok
  double precision anL, anR, ttL, ttR, SL, SR
  double precision uL_star, uR_star, umL_star, umR_star, u_star, p_star, pu_star
  double precision consLL(5),consRR(5),starLL(5),starRR(5), FL(5), FR(5),FF(5)

  double precision eps1

  if(keL.eq.keR) then
     print *, " + in RiemBN_flux_ "
     print *, " + forbiden case: keL, keR=", keL, keR
     pause
     stop
  end if

  eos_type = int(eos_stat(31,keL)+0.1)
  if(eos_type.eq.1) then
     GamL = eos_stat(1,keL)
  else
     print *, " + in RiemBN_flux, "
     print *, " + eos_type should be 1, but is ", eos_type
     pause
     stop
  end if

  eos_type = int(eos_stat(31,keR)+0.1)
  if(eos_type.eq.1) then
     GamR = eos_stat(1,keR)
  else
     print *, " + in RiemBN_flux, "
     print *, " + eos_type should be 1, but is ", eos_type
     pause
     stop
  end if

!  call eos_press(dL,eL,pL,cL,ttL,tmp,keL)
!  call eos_press(dR,eR,pR,cR,ttR,tmp,keR)

  flux_la(:) = 0.d0
  flux_eu(:) = 0.d0

  SL = dmin1(vxL-cL,vxR-cR)
  SR = dmax1(vxL+cL,vxR+cR)

  anL = dL*(vxL-SL)
  anR = dR*(SR-vxR)

  u_star = (pL-pR+anL*vxl+anR*vxR)/(anL+anR)
  p_star = pR+anR*(u_star-vxR)

  if(u_star.lt.SL.or.u_star.gt.SR) then
     print *, " + in RiemBN_flux: "
     print *, " + u_star, SL, SR= ",u_star, SL, SR
     pause
!     stop
  end if

  
  consLL(1) = dL
  consLL(2) = dL*vxL
  consLL(3) = dL*vyL
  consLL(4) = dL*vzL
  tmp = (vxL**2.d0+vyL**2.d0+vzL**2.d0)/2.d0
  consLL(5) = dL*(eL+tmp)  

  consRR(1) = dR
  consRR(2) = dR*vxR
  consRR(3) = dR*vyR
  consRR(4) = dR*vzR
  tmp = (vxR**2.d0+vyR**2.d0+vzR**2.d0)/2.d0
  consRR(5) = dR*(eR+tmp)

  starLL(1)=consLL(1)*(SL-vxL)/(SL-u_star)
  starLL(2)=starLL(1)*u_star
  starLL(3)=starLL(1)*vyL
  starLL(4)=starLL(1)*vzL
  tmp = (u_star-vxL)*(u_star+pL/dL/(SL-vxL))
  starLL(5)=starLL(1)*(consLL(5)/dL+tmp)

  starRR(1)=consRR(1)*(SR-vxR)/(SR-u_star)
  starRR(2)=starRR(1)*u_star
  starRR(3)=starRR(1)*vyR
  starRR(4)=starRR(1)*vzR
  tmp = (u_star-vxR)*(u_star+pR/dR/(SR-vxR))
  starRR(5)=starRR(1)*(consRR(5)/dR+tmp)

  
  ! calculate flux_la
  flux_la(0) = u_star
  flux_la(1) = 0.d0
  flux_la(2) = p_star
  flux_la(3) = 0.d0
  flux_la(4) = 0.d0
  flux_la(5) = p_star*u_star

  ! calculate flux_eu
  if(u_star.gt.0.d0) then   

     flux_eu(0)= u_star

     if(SL.le.0.d0) then
        flux_eu(1:5)= u_star*starLL(1:5)
     else
        flux_eu(1:5)= (u_star-SL)*starLL(1:5)+SL*consLL(1:5)
     endif

     flux_eu(6)= 1.d0     

  else 

     flux_eu(0)= u_star

     if(SR.ge.0.d0) then
        flux_eu(1:5)= u_star*starRR(1:5)
     else
        flux_eu(1:5)= (u_star-SR)*starRR(1:5)+SR*consRR(1:5)
     endif

     flux_eu(6)= -1.d0          

  end if
  
  
end subroutine RiemBN_flux_u2ph
