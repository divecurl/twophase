subroutine HLL_flux_u2ph(flux,u_movef,dL,eL,pL,cL,vxL,vyL,vzL,keL,dR,eR,pR,cR,vxR,vyR,vzR,keR)

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
  double precision consLL(5),consRR(5),cons_hll(5),FL(5),FR(5),FF(5)

!  call eos_press(dL,eL,pL,cL,ttL,tmp,keL)
!  call eos_press(dR,eR,pR,cR,ttR,tmp,keR)

  if(keL.ne.keR) then
     print *, " + in interface_slide_flux"
     print *, " + forbiden case: keL, keR=", keL, keR
     pause
     stop
  end if

  SL = dmin1(vxL-cL,vxR-cR)
  SR = dmax1(vxL+cL,vxR+cR)

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


  FL(1:5) = consLL(1:5)*vxL
  FL(2) = FL(2)+pL
  FL(5) = FL(5)+pL*vxL

  FR(1:5) = consRR(1:5)*vxR
  FR(2) = FR(2)+pR
  FR(5) = FR(5)+pR*vxR

  cons_hll(1:5)= (SR*consRR(1:5)-SL*consLL(1:5)+ FL(1:5) -FR(1:5))/(SR-SL)

  if(SL-u_movef.ge.0.d0) then
     FF(1:5)= FL(1:5)- consLL(1:5)*u_movef
  else if(SR-u_movef.le.0.d0) then
     FF(1:5)= FR(1:5)- consRR(1:5)*u_movef
  else 
     FF(1:5)=(SR*FL(1:5) -SL*FR(1:5)+ SL*SR*(consRR(1:5)-consLL(1:5)))/(SR-SL)
     FF(1:5)= FF(1:5)- cons_hll(1:5)*u_movef
  end if

  flux(1:5) = FF(1:5)

end subroutine HLL_flux_u2ph

  

