subroutine calcu_del_vx_nod_u2ph(nfirst,nlast,ngc)

  use compute_method_choosing,only: Nmat
  use Local_variables

  implicit none

  integer :: nfirst, nlast, ngc

  integer :: i,j,k,m,n, kL, kR


  double precision tmp, tmp1, tmp2

  double precision slopL, slopC, slopR

  do i=nfirst-ngc+1,nlast+ngc-2

     slopL = (vx_1d(2,i) - vx_1d(2,i-1))*2.d0/(dx_1d(i)+dx_1d(i-1))
     slopC = (vx_1d(2,i+1) - vx_1d(2,i))*2.d0/(dx_1d(i+1)+dx_1d(i))
     slopR = (vx_1d(2,i+2) - vx_1d(2,i+1))*2.d0/(dx_1d(i+2)+dx_1d(i+1))     

     if(slopC*slopL.le.0.d0.or.slopC*slopR.le.0.d0) then

        del_vx_nodL(i) = 0.d0
        del_vx_nodR(i) = 0.d0

     else  

        tmp = dmin1(dabs(slopL),dmin1(dabs(slopC),dabs(slopR)) )

        if(slopC.lt.0.d0) tmp = - tmp

        del_vx_nodL(i) = tmp * dx_1d(i)/2.d0

        del_vx_nodR(i) = tmp * dx_1d(i+1)/2.d0

     end if

  end do

end subroutine calcu_del_vx_nod_u2ph
