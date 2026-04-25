subroutine HLLC_flux_u2ph(flux,u_movef,dL,eL,pL,cL,vxL,vyL,vzL,keL,dR,eR,pR,cR,vxR,vyR,vzR,keR)

  use compute_method_choosing,only: flux_solver,eos_stat
  implicit none

  double precision flux(5),u_movef

  double precision dL,eL,pL,cL,vxL,vyL,vzL,dR,eR,pR,cR,vxR,vyR,vzR

  double precision av_L(2), av_R(2)

  integer :: keL, keR

  double precision tmp, tmp1, tmp2
  integer :: itmp
  double precision anL, anR, ttL, ttR, thetaL, thetaR, SL, SR
  double precision u_star, p_star, v_star, Sxy_star
  double precision consLL(5),consRR(5),starLL(5),starRR(5), FL(5), FR(5),FF(5)
  ! flux_a(1): u_star
  ! flux_a(2): p_star
  ! flux_a(3): u_star*p_star
  !    ! flux_a(4): v_star
  !    ! flux_a(5): Sxy_star
  !    ! flux_a(6): v_star*Sxy_star
  !    !!! flux_a(7:11) for sliding material
  ! flux_a(7): mas flux
  ! flux_a(8): xmom flux
  ! flux_a(9): ymom flux
  ! flux_a(10): zmom flux
  ! flux_a(11): energy flux
  !    ! flux_a(12): v_starL ! 用于滑移界面
  !    ! flux_a(13): v_starR ! 用于滑移界面

!  call eos_press(dL,eL,pL,cL,ttL,tmp,keL)
!  call eos_press(dR,eR,pR,cR,ttR,tmp,keR)

  if(keL.ne.keR) then
     print *, " + in interface_slide_flux"
     print *, " + forbiden case: keL, keR=", keL, keR
     pause
     stop
  end if

  flux(:) = 0.d0

  SL = dmin1(vxL-cL,vxR-cR)
  SR = dmax1(vxL+cL,vxR+cR)

  anL = dL*(vxL-SL)
  anR = dR*(SR-vxR)

  u_star = (pL-pR+anL*vxl+anR*vxR)/(anL+anR)
  p_star = pR+anR*(u_star-vxR)

  if(u_star.lt.SL.or.u_star.gt.SR) then
     print *, " + in slide_flux: "
     print *, " + u_star, SL, SR= ",u_star, SL, SR
     pause
     stop
  end if

  flux(1)= u_star
  flux(2)= p_star
  flux(3)= p_star*u_star

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

  FL(1:5) = consLL(1:5)*vxL
  FL(2) = FL(2)+pL
  FL(5) = FL(5)+pL*vxL

  FR(1:5) = consRR(1:5)*vxR
  FR(2) = FR(2)+pR
  FR(5) = FR(5)+pR*vxR

  if(SL-u_movef.ge.0.d0) then
     FF(1:5)= FL(1:5)- consLL(1:5)*u_movef
  else if(SR-u_movef.le.0.d0) then
     FF(1:5)= FR(1:5)- consRR(1:5)*u_movef
  else if(u_star-u_movef.ge.0.d0) then
     FF(1:5)= FL(1:5)+SL*(starLL(1:5)-consLL(1:5))-starLL(1:5)*u_movef
  else
     FF(1:5)= FR(1:5)+SR*(starRR(1:5)-consRR(1:5))-starRR(1:5)*u_movef
  end if

  flux(1:5) = FF(1:5)

end subroutine HLLC_flux_u2ph
  

