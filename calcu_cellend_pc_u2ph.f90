subroutine calcu_cellend_pc_u2ph(nfirst,nlast,ngc)

  use compute_method_choosing,only: Nmat
  use Local_variables

  implicit none

  integer :: nfirst, nlast, ngc

  integer :: i,j,k,m,n, kL, kR


  double precision tmp, tt_loc


  do i=nfirst-ngc,nlast+ngc-1

     do m=1,2
        call eos_press(dR_1d(m,i),eR_1d(m,i),pR_1d(m,i),cR_1d(m,i),tt_loc,tmp,ke_1d(m,i))
        call eos_press(dL_1d(m,i+1),eL_1d(m,i+1),pL_1d(m,i+1),cL_1d(m,i+1),tt_loc,tmp,ke_1d(m,i+1))
     end do

  end do

end subroutine calcu_cellend_pc_u2ph
