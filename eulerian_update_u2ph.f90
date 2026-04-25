subroutine eulerian_update_u2ph(dt, nfirst,nlast,ngc)

  use compute_method_choosing,only: Nmat,eos_stat
  use Local_variables

  implicit none

  double precision dt
  integer :: nfirst,nlast,ngc

  integer :: i,j,k,m,n

  double precision vol0,vol1,vol_mtr(5),vof_loc(5),mas_loc(5),xmo_loc(5),ymo_loc(5),zmo_loc(5),eng_loc(5)

  double precision vx_loc(5),vy_loc(5),vz_loc(5),dd_bck(5),ee_bck(5),pp_bck(5),eng_bck(5)

!  double precision  D_vofL(2),D_masL(2),D_xmoL(2),D_ymoL(2),D_zmoL(2),D_engL(2)

!  double precision  D_vofR(2),D_masR(2),D_xmoR(2),D_ymoR(2),D_zmoR(2),D_engR(2)

!  double precision  D_vofM(2),D_masM(2),D_xmoM(2),D_ymoM(2),D_zmoM(2),D_engM(2)

  double precision tmp, tmp1, tmp2, tmp3, ux_x, uy_x, anL, anR

  double precision weig(1:2)

  double precision u_act, p_act, pu_act


  do i=nfirst-1, nlast

     tmp1 = Mu_nod(i)
     tmp2= 1.d0-Mu_nod(i)

     do m=1,2
        fluxR(1:6,m,i) = tmp1*fluxR(1:6,m,i)+tmp2*fluxR2(1:6,m,i)
        fluxL(1:6,m,i+1) = tmp1*fluxL(1:6,m,i+1)+tmp2*fluxL2(1:6,m,i+1)
     enddo

  enddo

  do i=nfirst, nlast

     tmp1 = Mu_cell(i)
     tmp2= 1.d0-Mu_cell(i)
     
     do m=1,2
        fluxC(1:6,m,i) = tmp1*fluxC(1:6,m,i)+tmp2*fluxC2(1:6,m,i)
     enddo

  enddo


  do i=nfirst, nlast

     vol0= dx_1d(i)
           
     do m=1,2     

        D_vof_1d(m,i) = (fluxL(1,m,i) + fluxR(1,m,i) + fluxC(1,m,i))*dt/vol0 
        D_mas_1d(m,i) = (fluxL(2,m,i) + fluxR(2,m,i))*dt/vol0
        D_xmo_1d(m,i) = (fluxL(3,m,i) + fluxR(3,m,i) + fluxC(3,m,i))*dt/vol0
        D_ymo_1d(m,i) = (fluxL(4,m,i) + fluxR(4,m,i))*dt/vol0 
        D_zmo_1d(m,i) = (fluxL(5,m,i) + fluxR(5,m,i))*dt/vol0 
        D_eng_1d(m,i) = (fluxL(6,m,i) + fluxR(6,m,i) + fluxC(6,m,i))*dt/vol0

     end do

  end do

end subroutine eulerian_update_u2ph
