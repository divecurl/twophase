subroutine calcu_weig_mu_u2ph(weig_mu ,vof)

  use compute_method_choosing,only: Vof_crit, Weig_mu_min

  implicit none

  double precision weig_mu, vof

  if(Vof_crit(1).ge.Vof_crit(2) .or.Vof_crit(1).le.0.d0 .or.Vof_crit(2).ge.1.d0 ) then
      print *, "+ inreasonable Vof_crit(1:2) =", Vof_crit(1:2)
      pause
      stop
  endif    
  
  if(vof.le.Vof_crit(1)) then
     weig_mu= 0.d0
  else if(vof.ge.Vof_crit(2)) then
     weig_mu= 1.d0
  else
     weig_mu= (vof- vof_crit(1))/(vof_crit(2)-vof_crit(1))
  end if

  weig_mu = dmax1(Weig_mu_min,weig_mu)

end subroutine calcu_weig_mu_u2ph
