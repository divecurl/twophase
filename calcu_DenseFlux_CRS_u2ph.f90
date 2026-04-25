subroutine calcu_DenseFlux_CRS_u2ph(nfirst,nlast,ngc)

  use compute_method_choosing,only: Nmat
  use Local_variables

  implicit none

  integer :: nfirst, nlast,ngc

  integer :: i,j,k,m,n, kL, kR

  double precision flux_loc(0:5) ! to be solved by traditional Riemann solver, e.g. HLLC
  ! flux_loc(0) : 体积通量
  ! flux_loc(1) : 质量通量
  ! flux_loc(2) : x-动量通量
  ! flux_loc(3) : y-动量通量
  ! flux_loc(4) : z-动量通量
  ! flux_loc(5) : 能量通量

  double precision flux_la(0:5) ! Lagrangian part of flux by two phase SyBN solver 
  ! flux_la(0) : 体积通量
  ! flux_la(1) : 质量通量
  ! flux_la(2) : x-动量通量
  ! flux_la(3) : y-动量通量
  ! flux_la(4) : z-动量通量
  ! flux_la(5) : 能量通量

  double precision flux_eu(0:6) ! convection part of flux by two phase SyBN solver

  double precision vof_gas_loc, vof_solid_loc, vof_int

  double precision theL, theR, tmp, weig(2), p_act, u_act, pu_act

  do i=nfirst-ngc,nlast+ngc-1

     if(cell_type(i).eq.1.and.cell_type(i+1).eq.1) then
        mif_type(i)=0
     else if(cell_type(i).eq.1.and.cell_type(i+1).eq.2) then
        mif_type(i)=1
     else if(cell_type(i).eq.2 .and. cell_type(i+1).eq.1) then
        mif_type(i)=2
     else if(cell_type(i).eq.2.and.cell_type(i+1).eq.2) then
        if(dabs(vofR_1d(2,i)- vofL_1d(2,i+1)).lt.1.d-11) then
           mif_type(i)=3 
        else if(vofR_1d(2,i).lt.vofL_1d(2,i+1)) then
           mif_type(i)=4
        else
           mif_type(i)=5 
        endif
     else
        print *, "+ in calcu_flux_sybn :"
        print *, "+ this case of mif_type is not included !"
        continue
        pause
        stop
     end if

  end do


  do i=nfirst-ngc+1,nlast+ngc-2

     fluxR(:,:,i) = 0.d0
     fluxL(:,:,i+1) = 0.d0

     ! gas - gas interaction flux
     flux_loc(0:5) = 0.d0
     kL= 1
     kR= 1
!     call HLLC_flux_sybn(flux_loc(1:5),0.d0,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))
     call HLL_flux_u2ph(flux_loc(1:5),0.d0,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))

     if(mif_type(i).eq.0) then
        vof_gas_loc = 1.d0
     else 
        vof_gas_loc = dmin1(vofR_1d(1,i), vofL_1d(1,i+1))
     end if

     fluxR(2:6,1,i)= - flux_loc(1:5)*vof_gas_loc
     fluxL(2:6,1,i+1)= flux_loc(1:5)*vof_gas_loc


     ! solid - gas interactions
     if(mif_type(i).eq.0) then

        !               |   |   |  case(0)
        !               |   |   |

        print *, " this case of mif_type(i) is not included !"
        print *, " +  i, mif_type(i) = ", i, mif_type(i)        
        pause
        stop

     else if(mif_type(i).eq.1) then

        !               |   |___|  case(1)
        !               |   |   |

        print *, " this case of mif_type(i) is not included !"
        print *, " +  i, mif_type(i) = ", i, mif_type(i)        
        pause
        stop

     else if(mif_type(i).eq.2) then

        !               |___|   |  case(2)
        !               |   |   |
        print *, " this case of mif_type(i) is not included !"
        print *, " +  i, mif_type(i) = ", i, mif_type(i)        
        pause
        stop

     else if(mif_type(i).eq.3) then

        !               |___|___|   
        !               |   |   |  case(3)

        ! do nothing

     else if(mif_type(i).eq.4) then

        !               |   |___|   
        !               |___|   |  case(4)
        !               |   |   |


        flux_la(0:5) = 0.d0
        flux_eu(0:6) = 0.d0

        kL= 1
        kR= 2

        theL = cR_1d(kL,i)*dR_1d(kL,i)+cL_1d(kL,i+1)*dL_1d(kL,i+1)
        theR = cR_1d(kR,i)*dR_1d(kR,i)+cL_1d(kR,i+1)*dL_1d(kR,i+1)
        tmp = theL+theR
        theL = theL/tmp
        theR = theR/tmp

        ! BN model
        theL= 0.d0
        theR= 1.d0
        
        call SyBN_flux_u2ph(flux_la(0:5),flux_eu(0:6),theL,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),theR,dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))

        vof_int = vofR_1d(1,i)- vofL_1d(1,i+1)

        fluxR(1,1,i)= fluxR(1,1,i) + flux_la(0)*vof_int
        fluxR(2:6,1,i)= fluxR(2:6,1,i) - flux_la(1:5)*vof_int 
        
        fluxL(1,2,i+1)= fluxL(1,2,i+1) - flux_la(0)*vof_int
        fluxL(2:6,2,i+1)= fluxL(2:6,2,i+1) + flux_la(1:5)*vof_int 
        
        if(dabs(flux_eu(6)).lt.0.9d0) then

           print *, "+ inreasonable flux_eu(6): =", flux_eu(6)
           pause
           stop

        else if(flux_eu(6).gt.0) then

           fluxR(1,1,i)= fluxR(1,1,i) - flux_eu(0)*vof_int
           fluxR(2:6,1,i)= fluxR(2:6,1,i) - flux_eu(1:5)*vof_int 
           
           fluxL(1,1,i+1)= fluxL(1,1,i+1) + flux_eu(0)*vof_int
           fluxL(2:6,1,i+1)= fluxL(2:6,1,i+1) + flux_eu(1:5)*vof_int 

        else        

           fluxR(1,2,i)= fluxR(1,2,i) - flux_eu(0)*vof_int  ! flux_eu(0) < 0
           fluxR(2:6,2,i)= fluxR(2:6,2,i) - flux_eu(1:5)*vof_int 
           
           fluxL(1,2,i+1)= fluxL(1,2,i+1) + flux_eu(0)*vof_int
           fluxL(2:6,2,i+1)= fluxL(2:6,2,i+1) + flux_eu(1:5)*vof_int 

        end if

     else if(mif_type(i).eq.5) then     

        !               |___|   |   
        !               |   |___|  case(5)
        !               |   |   |

        flux_la(0:5) = 0.d0
        flux_eu(0:6) = 0.d0

        kL= 2
        kR= 1

        theL = cR_1d(kL,i)*dR_1d(kL,i)+cL_1d(kL,i+1)*dL_1d(kL,i+1)
        theR = cR_1d(kR,i)*dR_1d(kR,i)+cL_1d(kR,i+1)*dL_1d(kR,i+1)
        tmp = theL+theR
        theL = theL/tmp
        theR = theR/tmp

        ! BN model
        theL = 1.d0
        theR = 0.d0
        
        call SyBN_flux_u2ph(flux_la(0:5),flux_eu(0:6),theL,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),theR,dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))

        vof_int = dabs(vofR_1d(1,i)- vofL_1d(1,i+1))

        fluxR(1,2,i)= fluxR(1,2,i) + flux_la(0)*vof_int
        fluxR(2:6,2,i)= fluxR(2:6,2,i) - flux_la(1:5)*vof_int 
        
        fluxL(1,1,i+1)= fluxL(1,1,i+1) - flux_la(0)*vof_int
        fluxL(2:6,1,i+1)= fluxL(2:6,1,i+1) + flux_la(1:5)*vof_int 
        
        if(dabs(flux_eu(6)).lt.0.9d0) then

           print *, "+ inreasonable flux_eu(6): =", flux_eu(6)
           pause
           stop

        else if(flux_eu(6).gt.0) then

           fluxR(1,2,i)= fluxR(1,2,i) - flux_eu(0)*vof_int
           fluxR(2:6,2,i)= fluxR(2:6,2,i) - flux_eu(1:5)*vof_int 
           
           fluxL(1,2,i+1)= fluxL(1,2,i+1) + flux_eu(0)*vof_int
           fluxL(2:6,2,i+1)= fluxL(2:6,2,i+1) + flux_eu(1:5)*vof_int 

        else        

           fluxR(1,1,i)= fluxR(1,1,i) - flux_eu(0)*vof_int  ! flux_eu(0) < 0
           fluxR(2:6,1,i)= fluxR(2:6,1,i) - flux_eu(1:5)*vof_int 
           
           fluxL(1,1,i+1)= fluxL(1,1,i+1) + flux_eu(0)*vof_int
           fluxL(2:6,1,i+1)= fluxL(2:6,1,i+1) + flux_eu(1:5)*vof_int 

        end if

     end if

     ! solid -solid interaction

     flux_loc(0:5) = 0.d0
     kL= 2
     kR= 2
!     call HLLC_flux_sybn(flux_loc(1:5),0.d0,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))
     call HLL_flux_u2ph(flux_loc(1:5),0.d0,dR_1d(kL,i),eR_1d(kL,i),pR_1d(kL,i),cR_1d(kL,i),vxR_1d(kL,i),vyR_1d(kL,i),vzR_1d(kL,i),ke_1d(kL,i),dL_1d(kR,i+1),eL_1d(kR,i+1),pL_1d(kR,i+1),cL_1d(kR,i+1),vxL_1d(kR,i+1),vyL_1d(kR,i+1),vzL_1d(kR,i+1),ke_1d(kR,i+1))

     vof_solid_loc = dmin1(vofR_1d(2,i), vofL_1d(2,i+1))

     fluxR(2:6,2,i)= fluxR(2:6,2,i)  - flux_loc(1:5)*vof_solid_loc
     fluxL(2:6,2,i+1)= fluxL(2:6,2,i+1) + flux_loc(1:5)*vof_solid_loc

  end do

  
  ! cell center
  do i=nfirst-1,nlast+1

     fluxC(:,:,i) = 0.d0
      
     weig(1) = cc_1d(1,i)*dd_1d(1,i)
     weig(2) = cc_1d(2,i)*dd_1d(2,i)

     tmp = weig(1) +weig(2)
     weig(:)= weig(:)/tmp

    ! BN model
     weig(1) = 0.d0
     weig(2) = 1.d0
     
     u_act = weig(1)*vx_1d(1,i)+weig(2)*vx_1d(2,i)
     p_act = weig(2)*pp_1d(1,i)+weig(1)*pp_1d(2,i)
     pu_act = weig(2)*pp_1d(1,i)*vx_1d(2,i)+weig(1)*pp_1d(2,i)*vx_1d(1,i)

     do m=1,2

        tmp = vofL_1d(m,i) -  vofR_1d(m,i)
        fluxC(1,m,i) = u_act*tmp
        fluxC(3,m,i) = - p_act*tmp
        fluxC(6,m,i) = - pu_act*tmp

     end do

  end do

end subroutine calcu_DenseFlux_CRS_u2ph
