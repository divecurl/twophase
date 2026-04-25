subroutine dilute_flux5_u2ph(fluxL,fluxR,dL,eL,pL,cL,vxL,vyL,vzL,vofL,ppL_gas_in,uuL_solid_in,del_uL_solid,keL,&
    dR,eR,pR,cR,vxR,vyR,vzR,vofR,ppR_gas_in,uuR_solid_in,del_uR_solid,keR)
  
        !               |___|   |   
        !               |   |___|  case(5)
        !               |   |   |

  use compute_method_choosing,only: Mu_art_vis
  
  implicit none

  double precision fluxL(10,2), fluxR(10,2)
  ! fluxL_loc(1): 边界处气固耦合对左侧网格体积变化贡献(包含fluxL(7)的贡献) 
  ! fluxL_loc(2): 边界处气固耦合对左侧网格质量变化贡献 
  ! fluxL_loc(3): 边界处气固耦合对左侧网格动量变化贡献 x-direction
  ! fluxL_loc(4): 边界处气固耦合对左侧网格动量变化贡献 y-direction
  ! fluxL_loc(5): 边界处气固耦合对左侧网格动量变化贡献 z-direction
  ! fluxL_loc(6): 边界处气固耦合对左侧网格总能变化贡献(包含fluxL(7)的贡献)
  ! fluxL_loc(7): 边界处气固耦合固体不可压性引起左侧网格气固体积交换量
  ! fluxL_loc(8): 左侧网格向右侧网格输运的体积份额
  ! fluxL_loc(9): 左侧网格向右侧网格输运的固体体积份额
  ! fluxR_loc(1:9): 右侧网格的上述相应量

  double precision dL(2),eL(2),pL(2),cL(2),vxL(2),vyL(2),vzL(2),vofL(2)

  double precision ppL_gas_in, uuL_solid_in, del_uL_solid

  double precision dR(2),eR(2),pR(2),cR(2),vxR(2),vyR(2),vzR(2),vofR(2)

  double precision ppR_gas_in, uuR_solid_in, del_uR_solid

  integer :: keL(2), keR(2) 

  double precision uuL_solid, uuR_solid

  
!local variables
  double precision UL(5,2), UR(5,2)

  double precision ppL_gas, ppR_gas, pp_gas_itf
  double precision ppL_gas_vis, ppR_gas_vis
  double precision vof_cc(3)

  double precision tmp, tmp1, tmp2, p_sum, limL, limR
  double precision anL, anR, thetaL, thetaR, pppL,pppR  

  double precision u_itf
  integer :: i,j,k,m

  uuL_solid = uuL_solid_in
  uuR_solid = uuR_solid_in  
  
  ppL_gas= ppL_gas_in
  ppR_gas= ppR_gas_in
  
  pp_gas_itf= pR(1)

  !  ppL_gas= ppL_gas_in
  !  ppR_gas= ppR_gas_in
  
  !------------------------------------------------------------
  
  vof_cc(1)= dmin1(vofL(1),vofR(1))
  vof_cc(2)= vofR(1)-vofL(1)

  if(vof_cc(2).le.-1.d11) then
     print *, "+ in _flux5 inreasonable vof_cc(2): =", vof_cc(2)
     pause
     stop
  endif
  
  vof_cc(3)= dmin1(vofL(2),vofR(2))

! ----------------- dilute twophase solver ---------------
  
  fluxL(:,:)= 0.d0
  fluxR(:,:)= 0.d0

     anL = vofL(2)
     anR = vofR(2)
     thetaL = anL/(anL+anR)
     thetaR = anR/(anL+anR)

     u_itf= thetaL*uuL_solid+thetaR*uuR_solid

!     if(uuL_solid.gt.u_itf) then       
!        limL = 1.d0
!        tmp1 = uuL_solid-u_itf
!        if(tmp1*del_uL_solid.lt.1.d-12) then
!           tmp =dmin1(1.d0, -del_uL_solid/tmp1)
!           limL = 1.d0-tmp
!        endif
!        anL = dL(2)*cL(2)
!        tmp = anL*(uuL_solid-u_itf)*Mu_art_vis*limL
!        ppL_gas_vis = ppL_gas + tmp
!     else
        ppL_gas_vis = ppL_gas 
!     end if

!     if(u_itf.gt.uuR_solid) then       
!        limR = 1.d0
!        tmp2 = u_itf - uuR_solid
!        if(tmp2*del_uR_solid.lt.1.d-12) then        
!           tmp =dmin1(1.d0, -del_uR_solid/tmp2)          
!           limR = 1.d0-tmp
!        endif

!        anR = dR(2)*cR(2)
!        tmp = anR*(u_itf - uuR_solid)*Mu_art_vis*limR
!        ppR_gas_vis = ppR_gas + tmp
!     else
        ppR_gas_vis = ppR_gas 
!     end if

     p_sum= pp_gas_itf*vof_cc(2)+ppR_gas_vis*vofR(2)-ppL_gas_vis*vofL(2)

! flux on left
     fluxL(7,1) = (u_itf-uuL_solid)*vofL(2)
     fluxL(1,1)= fluxL(7,1)
     fluxL(3,1) = 0.d0
     fluxL(6,1) = -ppL_gas*fluxL(7,1) 

     fluxL(7,2) = -fluxL(7,1)
     fluxL(1,2)= u_itf*vofL(2)+fluxL(7,2)
     fluxL(3,2)= -ppL_gas_vis*vofL(2) - thetaL*p_sum
     tmp = -ppL_gas_vis*vofL(2)*u_itf-thetaL*p_sum*uuL_solid 
     fluxL(6,2)= tmp- ppL_gas*fluxL(7,2)

! flux on the right
     fluxR(7,1) = -(u_itf-uuR_solid)*vofR(2)
     fluxR(1,1) = - u_itf*vof_cc(2)+fluxR(7,1)
     fluxR(3,1)= pp_gas_itf*vof_cc(2)
     fluxR(6,1)= pp_gas_itf*u_itf*vof_cc(2)-ppR_gas*fluxR(7,1)

     fluxR(7,2) = -fluxR(7,1)
     fluxR(1,2)= -u_itf*vofR(2)+fluxR(7,2)
     fluxR(3,2)= ppR_gas_vis*vofR(2)- thetaR*p_sum
     tmp = ppR_gas_vis*vofR(2)*u_itf-thetaR*p_sum*uuR_solid 
     fluxR(6,2)= tmp- ppR_gas*fluxR(7,2)
     
     ! artificial viscosity for velocity equilibium
     uuL_solid = vxL(2)
     uuR_solid = vxR(2)
          
!     if(uuL_solid.gt.uuR_solid) then

        anL = vofL(2)*dL(2)*cL(2)
        anR = vofR(2)*dR(2)*cR(2)

        tmp= (anL*uuL_solid+anR*uuR_solid)/(anL+anR)
        tmp1 = 1.d0
                
        fluxL(3,2) = fluxL(3,2) + anL*(tmp-uuL_solid)*tmp1*Mu_art_vis
        fluxL(6,2) = fluxL(6,2) + anL*(tmp-uuL_solid)*tmp1*tmp*Mu_art_vis

        fluxR(3,2) = fluxR(3,2) + anR*(tmp-uuR_solid)*tmp1*Mu_art_vis
        fluxR(6,2) = fluxR(6,2) + anR*(tmp-uuR_solid)*tmp1*tmp*Mu_art_vis
        
 !    end if
     
! convection part

     if(u_itf.gt.0.d0) then

        tmp1 = 0.d0
        tmp2 = u_itf*vofL(2)

        do m=1,2
           UL(1,m) = dL(m)
           UL(2,m) = UL(1,m)*vxL(m)
           UL(3,m) = UL(1,m)*vyL(m)
           UL(4,m) = UL(1,m)*vzL(m)
           tmp= vxL(m)**2.d0+vyL(m)**2.d0+vzL(m)**2.d0
           UL(5,m) = UL(1,m)*(eL(m)+tmp/2.d0)
        end do

!        fluxL(1,1) = fluxL(1,1) - tmp1
!        fluxL(2:6,1) = fluxL(2:6,1) - tmp1*UL(1:5,1)

        fluxL(1,2) = fluxL(1,2) - tmp2
        fluxL(2:6,2) = fluxL(2:6,2) - tmp2*UL(1:5,2)

!        fluxR(1,1) = fluxR(1,1) + tmp1
!        fluxR(2:6,1) = fluxR(2:6,1) + tmp1*UL(1:5,1)

        fluxR(1,2) = fluxR(1,2) + tmp2
        fluxR(2:6,2) = fluxR(2:6,2) + tmp2*UL(1:5,2)

     else

        tmp1 = u_itf*vof_cc(2)
        tmp2 = u_itf*vofR(2)

        do m=1,2
           UR(1,m) = dR(m)
           UR(2,m) = UR(1,m)*vxR(m)
           UR(3,m) = UR(1,m)*vyR(m)
           UR(4,m) = UR(1,m)*vzR(m)
           tmp= vxR(m)**2.d0+vyR(m)**2.d0+vzR(m)**2.d0
           UR(5,m) = UR(1,m)*(eR(m)+tmp/2.d0)
        end do

        fluxL(1,1) = fluxL(1,1) - tmp1
        fluxL(2:6,1) = fluxL(2:6,1) - tmp1*UR(1:5,1)

        fluxL(1,2) = fluxL(1,2) - tmp2
        fluxL(2:6,2) = fluxL(2:6,2) - tmp2*UR(1:5,2)

        fluxR(1,1) = fluxR(1,1) + tmp1
        fluxR(2:6,1) = fluxR(2:6,1) + tmp1*UR(1:5,1)

        fluxR(1,2) = fluxR(1,2) + tmp2
        fluxR(2:6,2) = fluxR(2:6,2) + tmp2*UR(1:5,2)

     end if


end subroutine dilute_flux5_u2ph


