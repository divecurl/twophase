subroutine one_dim_deltaU_u2ph(dt,nfirst,nlast,ngc,id_direction)
    
  use compute_method_choosing,only:Nmat,eos_stat,vof_cut_eps,Twophase_solver
  use Local_variables

  implicit none

  double precision :: dt
  integer :: nfirst, nlast, ngc, id_direction    
  integer :: i,j,k,k1,m,n,itmp,itmp1,itmp2, itmp3, index_ok
  
  double precision vol0, vol1, vol2, vtmp(10), tmp, tmp1, tmp2

  double precision vof_loc(2), dd_loc(2), ee_loc(2), pp_loc(2), cc_loc(2), tt_loc(2)
  
  do i=nfirst-ngc,nlast+ngc
     do m=1,Nmat
        if(ke_1d(m,i).ge.0) then
           cell_type(i) = m
        end if
     end do
  end do

  ! high order reconstruction
  call reconstruct_vof_sybn(vofL_1d(1:Nmat,nfirst-ngc:nlast+ngc),vofR_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  call reconstruct_sybn(dL_1d(1:Nmat,nfirst-ngc:nlast+ngc),dR_1d(1:Nmat,nfirst-ngc:nlast+ngc),dd_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  call reconstruct_sybn(eL_1d(1:Nmat,nfirst-ngc:nlast+ngc),eR_1d(1:Nmat,nfirst-ngc:nlast+ngc),ee_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  call reconstruct_sybn(vxL_1d(1:Nmat,nfirst-ngc:nlast+ngc),vxR_1d(1:Nmat,nfirst-ngc:nlast+ngc),vx_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  call reconstruct_sybn(vyL_1d(1:Nmat,nfirst-ngc:nlast+ngc),vyR_1d(1:Nmat,nfirst-ngc:nlast+ngc),vy_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  call reconstruct_sybn(vzL_1d(1:Nmat,nfirst-ngc:nlast+ngc),vzR_1d(1:Nmat,nfirst-ngc:nlast+ngc),vz_1d(1:Nmat,nfirst-ngc:nlast+ngc),vof_1d(1:Nmat,nfirst-ngc:nlast+ngc),dx_1d(nfirst-ngc:nlast+ngc),cell_type(nfirst-ngc:nlast+ngc),nfirst,nlast,ngc,Nmat)

  if(twophase_solver.eq.1) then  

     call calcu_cellend_pc_u2ph(nfirst,nlast,ngc)

     call calcu_Mu_u2ph(nfirst,nlast,ngc)

     call calcu_del_vx_nod_u2ph(nfirst,nlast,ngc)

     call calcu_DenseFlux_CRS_u2ph(nfirst,nlast,ngc)

     call calcu_DiluteFlux_u2ph(nfirst,nlast,ngc)

     call eulerian_update_u2ph(dt,nfirst,nlast,ngc)
      
  else if(twophase_solver.eq.2) then 

     call calcu_cellend_pc_u2ph(nfirst,nlast,ngc)

     call calcu_DenseFlux_RiemModel_u2ph(nfirst,nlast,ngc)

     call eulerian_update_Riem_u2ph(dt,nfirst,nlast,ngc)
     
  else

     print *, "+ this solver is not expected !"
     print *, "+ twophase_solver = ", twophase_solver
     pause
     stop

  endif

end subroutine one_dim_deltaU_u2ph

