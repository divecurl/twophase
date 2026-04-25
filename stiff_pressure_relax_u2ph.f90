subroutine stiff_pressure_relax_u2ph(vof,dd,in_e,pp,cc,tt,ke)

  use compute_method_choosing,only: eos_stat, vof_cut_eps

  implicit none

  double precision, dimension(2):: vof, dd, in_e, pp, cc, tt
  double precision, dimension(2):: vof_in, dd_in, pp_in, in_e_in, cc_in, tt_in
  integer :: ke(2)

  integer :: i,j,k,m,n

  double precision, dimension(2):: Gam, Eta, Pi, Cv

  double precision mas(2), ant(2), Pint, A, B, C(2), D

  double precision tmp

  vof_in = vof
  dd_in = dd
  in_e_in = in_e
  
  
  if(int(eos_stat(31,ke(1))+0.1).ne.1 .or. int(eos_stat(31,ke(2))+0.1).ne.1) then
     print *, "+ in stiff_pressure_relax: "
     print *, "+ only eos type 1 is permitted here"
     pause
     stop
  end if

  Gam(1) = eos_stat(1,ke(1))
  Pi(1) = eos_stat(2,ke(1))
  Eta(1) = eos_stat(3,ke(1))
!  Cv(1) = eos_stat(4,ke(1))

  Gam(2) = eos_stat(1,ke(2))
  Pi(2) = eos_stat(2,ke(2))
  Eta(2) = eos_stat(3,ke(2))
!  Cv(2) = eos_stat(4,ke(2))

  call eos_press(dd_in(1),in_e_in(1),pp_in(1),cc_in(1),tt_in(1),Cv(1),ke(1))
  call eos_press(dd_in(2),in_e_in(2),pp_in(2),cc_in(2),tt_in(2),Cv(2),ke(2))

  ant(1)= dd_in(1)*cc_in(1)
  ant(2)= dd_in(2)*cc_in(2)

  Pint= (ant(1)*pp_in(2)+ant(2)*pp_in(1))/(ant(1)+ant(2))

  A=1.d0+Gam(1)*vof_in(2)+Gam(2)*vof_in(1)
  C(1)=2.d0*Gam(1)*Pi(1)+(Gam(1)-1.d0)*Pint
  C(2)=2.d0*Gam(2)*Pi(2)+(Gam(2)-1.d0)*Pint
  B=C(1)*vof_in(2)+C(2)*vof_in(1)-(1.d0+Gam(2))*vof_in(1)*pp_in(1)-(1.d0+Gam(1))*vof_in(2)*pp_in(2) 
  D= -C(2)*vof_in(1)*pp_in(1) -C(1)*vof_in(2)*pp_in(2)

  tmp = B**2.d0-4.d0*A*D
  if(tmp.lt.0.d0) then
     print *, "+ in stiff_pressure_relax:"
     print *, "+ B*B-4*A*D <0 := ", B**2.d0-4.d0*A*D
     pause
     stop
  else
     pp= (-B+dsqrt(tmp))/2.d0/A
  end if

  vof(1) = (Gam(1)-1.d0)*pp(1)+2.d0*pp_in(1)+C(1)
  vof(1)=vof(1)*vof_in(1)/((Gam(1)+1.d0)*pp(1)+C(1))

  vof(2)= 1.d0- vof(1)

  if(dmin1(vof(1),vof(2)).lt. vof_cut_eps/10.d0.or.pp(1).lt.0.d0 ) then
     print *, "+ in stiff_pressure_relax:"
     print *, "+ vof(:) =", vof(:)
     print *, "+ pp(:) =", pp(:)
     pause
     stop
  end if

  mas(1:2)= dd_in(1:2)*vof_in(1:2)
  dd(1:2)= mas(1:2)/vof(1:2)
  
  in_e(1)= (pp(1)+Gam(1)*Pi(1))/(Gam(1)-1.d0)/dd(1) + Eta(1)
  in_e(2)= (pp(2)+Gam(2)*Pi(2))/(Gam(2)-1.d0)/dd(2) + Eta(2)

  call eos_press(dd(1),in_e(1),pp(1),cc(1),tt(1),Cv(1),ke(1))
  call eos_press(dd(2),in_e(2),pp(2),cc(2),tt(2),Cv(2),ke(2))

end subroutine stiff_pressure_relax_u2ph
