subroutine SyBN_flux_u2ph(flux_la,flux_eu,theL,dL,eL,pL,cL,vxL,vyL,vzL,keL,theR,dR,eR,pR,cR,vxR,vyR,vzR,keR)

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
  double precision uL_star, uR_star, umL_star, umR_star, u_int_star, p_int, pu_int_star
  double precision consLL(5),consRR(5),starLL(5),starRR(5), FL(5), FR(5),FF(5)

  double precision eps1

  if(keL.eq.keR) then
     print *, " + in SyBN_flux_ "
     print *, " + forbiden case: keL, keR=", keL, keR
     pause
     stop
  end if

  eos_type = int(eos_stat(31,keL)+0.1)
  if(eos_type.eq.1) then
     GamL = eos_stat(1,keL)
  else
     print *, " + in SyBN_flux, "
     print *, " + eos_type should be 1, but is ", eos_type
     pause
     stop
  end if

  eos_type = int(eos_stat(31,keR)+0.1)
  if(eos_type.eq.1) then
     GamR = eos_stat(1,keR)
  else
     print *, " + in SyBN_flux, "
     print *, " + eos_type should be 1, but is ", eos_type
     pause
     stop
  end if

!  call eos_press(dL,eL,pL,cL,ttL,tmp,keL)
!  call eos_press(dR,eR,pR,cR,ttR,tmp,keR)

  flux_la(:) = 0.d0
  flux_eu(:) = 0.d0

  p_int = theL*pR+theR*pL
  u_int_star = theL*vxL+theR*vxR

  SL = dmin1(vxL-cL,vxR-cR,u_int_star-cL)
  SR = dmax1(vxL+cL,vxR+cR,u_int_star+cR)
  

  anL = dL*(vxL-SL)
  anR = dR*(SR-vxR)
     
  uL_star = vxL+ (pL-p_int)/anL
  uR_star = vxR+ (p_int-pR)/anR

  pu_int_star= theR*pL*vxR+theL*pR*vxL  
  
  consLL(1) = dL
  consLL(2) = dL*vxL
  consLL(3) = dL*vyL
  consLL(4) = dL*vzL
  tmp = (vxL**2.d0+vyL**2.d0+vzL**2.d0)/2.d0
  consLL(5) = dL*(eL+tmp)  

  starLL(1) = consLL(1)*(SL-vxL)/(SL-u_int_star)
  starLL(2) = starLL(1)*uL_star
  starLL(3) = starLL(1)*vyL
  starLL(4) = starLL(1)*vzL
  starLL(5) = (consLL(5)*(SL-vxL)+(pu_int_star-pL*vxL))/(SL-u_int_star)

  consRR(1) = dR
  consRR(2) = dR*vxR
  consRR(3) = dR*vyR
  consRR(4) = dR*vzR
  tmp = (vxR**2.d0+vyR**2.d0+vzR**2.d0)/2.d0
  consRR(5) = dR*(eR+tmp)

  starRR(1) = consRR(1)*(SR-vxR)/(SR-u_int_star)
  starRR(2) = starRR(1)*uR_star
  starRR(3) = starRR(1)*vyR
  starRR(4) = starRR(1)*vzR
  starRR(5) = (consRR(5)*(SR-vxR)+(pu_int_star-pR*vxR))/(SR-u_int_star)


  u_int_star = theL*vxL+theR*vxR
  p_int = theL*pR+theR*pL  
  pu_int_star = theL*vxL*pR+theR*vxR*pL  
  
  ! calculate flux_la
  flux_la(0) = u_int_star
  flux_la(1) = 0.d0
  flux_la(2) = p_int
  flux_la(3) = 0.d0
  flux_la(4) = 0.d0
  flux_la(5) = pu_int_star

  ! calculate flux_eu
  if(u_int_star.gt.0.d0) then   
     flux_eu(0)= u_int_star
     flux_eu(1:5)= u_int_star*consLL(1:5)
     flux_eu(6)= 1.d0     
  else 
     flux_eu(0)= u_int_star
     flux_eu(1:5)= u_int_star*consRR(1:5)
     flux_eu(6)= -1.d0          
  end if
  
  
end subroutine SyBN_flux_u2ph
