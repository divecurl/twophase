subroutine calcu_Mu_u2ph(nfirst,nlast,ngc)

  use compute_method_choosing,only: Nmat
  use Local_variables

  implicit none

  integer :: nfirst, nlast, ngc

  integer :: i,j,k,m,n, kL, kR

  double precision tmp, tmp1, tmp2

  double precision slopL, slopC, slopR

  do i=nfirst-ngc,nlast+ngc-1

     tmp = dmax1(vofR_1d(2,i),vofL_1d(2,i+1))

     call calcu_weig_mu_u2ph(tmp1,tmp)

     Mu_nod(i) = tmp1

  end do

  do i=nfirst-ngc,nlast+ngc

     tmp = vof_1d(2,i)

     call calcu_weig_mu_u2ph(tmp1,tmp)

     Mu_cell(i) = tmp1

  end do
  
end subroutine calcu_Mu_u2ph
