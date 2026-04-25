subroutine calcu_DiluteFlux_u2ph(nfirst,nlast,ngc)

! calculate fluxL2, fluxR2, fluxC3

  use compute_method_choosing,only: Nmat
  use Local_variables

  implicit none

  integer :: nfirst, nlast,ngc

  integer :: i,j,k,m,n, kL, kR

  double precision flux_gas_loc(5) ! 气流通道中的气体flux

  double precision fluxL_loc(10,2), fluxR_loc(10,2) 
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

  double precision ucL(5,2), ucR(5,2)

  double precision vof_gas_loc, tmp

  double precision weig(1:2)

  double precision u_act, p_act, pu_act

  integer :: nfi_ac, nla_ac, itmp

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  nfi_ac= nfirst    
  nla_ac= nlast    

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    

  do i=nfirst-ngc,nlast+ngc-1

     if(vofR_1d(2,i).le.vofL_1d(2,i+1)) then
        mif_type(i)=4
     else
        mif_type(i)=5 
     endif

  end do


  do i=nfi_ac-2, nla_ac+2

     fluxL2(:,:,i) = 0.d0

     fluxR2(:,:,i) = 0.d0

     fluxC2(:,:,i) = 0.d0

  end do

  do i=nfi_ac-1, nla_ac

     ! gas - gas interaction flux
     flux_gas_loc(1:5) = 0.d0
     kL= 1
     kR= 1
     call HLL_flux_u2ph(flux_gas_loc(1:5),0.d0,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),&
     ke_1d(kL,i),dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))
     if(mif_type(i).eq.0) then
        vof_gas_loc = 1.d0
     else 
        vof_gas_loc = dmin1(vofR_1d(1,i), vofL_1d(1,i+1))
     end if
     flux_gas_2ph(:,i) =  flux_gas_loc(:)*vof_gas_loc


     fluxR2(1,1,i) = 0.d0
     fluxR2(2:6,1,i) = - flux_gas_2ph(1:5,i)

     fluxL2(1,1,i+1) = 0.d0
     fluxL2(2:6,1,i+1) = flux_gas_2ph(1:5,i)


     ! gas solid interaction

     if (Mu_nod(i).lt.1.d0-1.d-12) then     
     
         if(mif_type(i).eq.4) then

            !               |   |___|   
            !               |___|   |  case(4)
            !               |   |   |

            call dilute_flux4_u2ph(fluxL_loc,fluxR_loc,dR_1d(1:2,i),eR_1d(1:2,i),pR_1d(1:2,i),cR_1d(1:2,i),vxR_1d(1:2,i),vyR_1d(1:2,i),vzR_1d(1:2,i),vofR_1d(1:2,i),pp_1d(1,i),vx_1d(2,i),del_vx_nodL(i), ke_1d(1:2,i),dL_1d(1:2,i+1),eL_1d(1:2,i+1),pL_1d(1:2,i+1),cL_1d(1:2,i+1),vxL_1d(1:2,i+1),vyL_1d(1:2,i+1),vzL_1d(1:2,i+1),vofL_1d(1:2,i+1),pp_1d(1,i+1),vx_1d(2,i+1), del_vx_nodR(i),ke_1d(1:2,i+1))

        
         else if(mif_type(i).eq.5) then     

            !               |___|   |   
            !               |   |___|  case(5)
            !               |   |   |

            call dilute_flux5_u2ph(fluxL_loc,fluxR_loc,dR_1d(1:2,i),eR_1d(1:2,i),pR_1d(1:2,i),cR_1d(1:2,i),vxR_1d(1:2,i),vyR_1d(1:2,i),vzR_1d(1:2,i),vofR_1d(1:2,i),pp_1d(1,i),vx_1d(2,i),del_vx_nodL(i),ke_1d(1:2,i),dL_1d(1:2,i+1),eL_1d(1:2,i+1),pL_1d(1:2,i+1),cL_1d(1:2,i+1),vxL_1d(1:2,i+1),vyL_1d(1:2,i+1),vzL_1d(1:2,i+1),vofL_1d(1:2,i+1),pp_1d(1,i+1),vx_1d(2,i+1),del_vx_nodR(i),ke_1d(1:2,i+1))

         else

            print *, "+ this case of mif_type is not included,"
            print *, "+ mif_type(i), i =", mif_type(i), i
            pause
            stop

         end if

     else
         
         fluxL_loc(:,:) = 0.d0
         fluxR_loc(:,:) = 0.d0
     
     endif
     
     fluxR2(1:6,1:2,i)= fluxR2(1:6,1:2,i) + fluxL_loc(1:6,1:2)
     fluxL2(1:6,1:2,i+1)= fluxL2(1:6,1:2,i+1) + fluxR_loc(1:6,1:2)

  enddo

  do i=nfi_ac, nla_ac

     weig(1) = 0.d0
     weig(2) = 1.d0

     u_act = weig(1)*vx_1d(1,i)+weig(2)*vx_1d(2,i)
     p_act = weig(2)*pp_1d(1,i)+weig(1)*pp_1d(2,i)
     pu_act = weig(2)*pp_1d(1,i)*vx_1d(2,i)+weig(1)*pp_1d(2,i)*vx_1d(1,i)

     do m=1,2
         tmp = vofL_1d(m,i) -  vofR_1d(m,i)
         fluxC2(1,m,i) = u_act*tmp
         fluxC2(3,m,i) = - p_act*tmp
         fluxC2(6,m,i) = - pu_act*tmp
     enddo
     
  enddo

end subroutine calcu_DiluteFlux_u2ph
