subroutine eos_press(dd,ee,pp,cc,tempr,cv,mat_id)

  use compute_method_choosing,only: eos_stat
    
  implicit none

  integer :: mat_id

  double precision dd !  密度
  double precision ee ! 内能
  double precision cc ! 声速
  double precision pp ! 压力
  double precision tempr ! 温度
  double precision cv  ! d_ee/d_tpr
  
  double precision SG_Gam, SG_Pi, SG_Eta, SG_Cv

  double precision MG_Gam, MG_Rho, MG_c, MG_s, MG_R, x_loc, p_ref, e_ref, ff, df 

  double precision e00, e01, e02, e03, e04
  
  double precision rho_min, p_min, ine_min, c_min 

  integer :: eos_type
  
  eos_type = int(eos_stat(31,mat_id)+0.1)

  select case(eos_type)

  case(0)    

!      print *, " + in eos: "
!      print *, " + eos_type = ", eos_type
!      print *, " + do nothing ..."
      pp = 1.d-6
      cc = 1.d-6
      tempr = 1.d-6
      cv = 1.d0
      dd= 1.d-6
      ee= 1.d-6
      
  case(1)

     SG_Gam = eos_stat(1,mat_id)
     SG_Pi = eos_stat(2,mat_id)
     SG_Eta = eos_stat(3,mat_id)
     SG_Cv = eos_stat(4,mat_id)
     
     rho_min = eos_stat(17,mat_id)
     p_min = eos_stat(18,mat_id)
     ine_min = eos_stat(19,mat_id)
     c_min = eos_stat(20,mat_id)     

     if(dd.lt.rho_min .or. ee.lt.ine_min) then
        !print *, " in  eos_press_sound: "
        !print *, " rho_min, dd = ", rho_min, dd
        !print *, " ine_min, ee = ", ine_min, ee
        !print *, " mat_id, eos_type = ", mat_id, eos_type
        !continue
!        pause
     end if

     dd = dmax1(dd,rho_min)     
     ee = dmax1(ee,ine_min)
     
     !if(ee.ge.5.d0) then
     !    continue
     !endif    
     
     pp = (SG_Gam-1.d0)*dd*ee - (SG_Gam-1.d0)*SG_Eta*dd - SG_Gam*SG_Pi

     cv = SG_Cv
     
     tempr= (pp+SG_Pi)/cv/dd/(SG_Gam-1.d0)
     
     tempr = max(1.d-8,tempr)
     
     cc = SG_Gam*(pp+SG_Pi)/dd

     if(pp.lt.p_min .or. cc.lt.c_min**2.d0 ) then
        !print *, " in  eos_press_sound: "
        !print *, " p_min, pp = ", p_min, pp
        !print *, " c_min^2, cc^2 = ", c_min**2.d0, cc
        !print *, " mat_id, eos_type = ", mat_id, eos_type
        !continue
        !pause
        pp=dmax1(pp,p_min)
        cc=dmax1(cc, c_min**2.d0 )
        
     end if

     cc= sqrt(cc)

  case(2)

     MG_Gam = eos_stat(1,mat_id)
     MG_Rho = eos_stat(2,mat_id)
     MG_c = eos_stat(3,mat_id)
     MG_s = eos_stat(4,mat_id)
     MG_R = eos_stat(5,mat_id)

     rho_min = eos_stat(17,mat_id)
     p_min = eos_stat(18,mat_id)
     ine_min = eos_stat(19,mat_id)
     c_min = eos_stat(20,mat_id)     
 

     if(dd.lt.rho_min .or. ee.lt.ine_min) then
        print *, " in  eos_press_sound: "
        print *, " rho_min, dd = ", rho_min, dd
        print *, " ine_min, ee = ", ine_min, ee
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if       
     
     
!     if(dd.lt.rho_min .or. ee.lt.ine_min) then
!        print *, " in  eos_press_sound: "
!        print *, " rho_min, dd = ", rho_min, dd
!        print *, " ine_min, ee = ", ine_min, ee
!        print *, " mat_id, eos_type = ", mat_id, eos_type
!        continue
!        pause
!     end if

!ccccccccccccccccccccccccccccccccccccccccccccccc 计算温度
     
     x_loc= 1.d0- MG_Rho/dd
     if(1.d0-MG_s*x_loc.le.0.3d0) then
         print *, " + in eos_press: 1.d0-MG_s*x_loc is small = ", 1.d0-MG_s*x_loc
         pause
     endif
     
     e00= -900.d0*MG_R
     e01= MG_Gam*e00
     e02= (MG_c**2.d0+MG_Gam**2.d0*e00)/2.d0
     e03= (4.d0*MG_s*MG_c**2.d0+MG_Gam**3.d0*e00)/6.d0
     e04= (18.d0*MG_s**2.d0*MG_c**2.d0-2.d0*MG_Gam*MG_s*MG_c**2.d0+ MG_Gam**4.d0*e00)/24.d0
          
     e_ref= e00+e01*x_loc+e02*x_loc**2.d0+e03*x_loc**3.d0+e04*x_loc**4.d0

     cv = 3.d0*MG_R
     tempr= (ee-e_ref)/(cv+1.d-99)
     if(tempr.le.280.d0) then
         !print *, "+ in eos_press: tempr is small =", tempr
         !print *, "dd= ", dd         
!         pause
     endif    
     tempr= max(198d0,tempr)
    
!cccccccccccccccccccccccccccccccccccccccccccccc  计算压力     
     
     if(x_loc.ge.0.d0) then     
        ff= MG_Rho*MG_c**2.d0*x_loc*(1.d0-MG_Gam*x_loc/2.d0)/(1.d0-MG_s*x_loc)**2.d0
        df= MG_Rho*MG_c**2.d0*(1.d0+(MG_s-MG_Gam)*x_loc)/(1.d0-MG_s*x_loc)**3.d0
        df= MG_Rho/dd/dd*df
     else
        ff= MG_c**2.d0*(dd-MG_Rho)
        df= MG_c**2.d0
     endif  
        
     pp= MG_Gam*MG_Rho*ee + ff     
     cc= MG_Gam*MG_Rho/dd/dd*pp+df


     if(dd.lt.rho_min .or. ee.lt.ine_min) then
        print *, " in  eos_press_sound: "
        print *, " rho_min, dd = ", rho_min, dd
        print *, " ine_min, ee = ", ine_min, ee
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if         
     
     if(pp.lt.p_min*1.0001d0 .or. cc.lt.c_min**2.d0 ) then
        !print *, " in  eos_press_sound: "
        !print *, " p_min, pp, dd = ", p_min, pp, dd
        !print *, " c_min^2, cc^2 = ", c_min**2.d0, cc
        !print *, " mat_id, eos_type = ", mat_id, eos_type
        !continue
        !pause
        
        pp=dmax1(pp,p_min*1.0001)
        cc=dmax1(cc, c_min**2.d0 )
        
     end if

     cc= sqrt(cc)
    
     
     if(cc.ge.MG_c*3.d0) then 
        print *, " in  eos_press_sound: "
        print *, " p_min, pp, dd = ", p_min, pp, dd
        print *, " cc is so large = ",  cc
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
         
        cc= dmin1(cc, MG_c*3.d0)
     endif   
        
        
  case default
        
     print *, " in eos_ienergy: this is case of eos_type is not included ."
     pause
     stop

  end select
  
end subroutine eos_press


subroutine eos_ienergy(dd,pp,ee,mat_id)

  use compute_method_choosing,only : eos_stat

  implicit none

  integer :: mat_id

  double precision dd, ee, pp

  double precision SG_Gam, SG_Pi, SG_Eta, SG_Cv

  double precision MG_Gam, MG_Rho, MG_c, MG_s, MG_R, x_loc, ff, df 

  double precision rho_min, p_min, ine_min, c_min 

  integer :: eos_type
  
  eos_type = int(eos_stat(31,mat_id) +0.1)

  select case(eos_type)

  case(0)
      ee =1.d-6
      
  case(1)

     SG_Gam = eos_stat(1,mat_id)
     SG_Pi = eos_stat(2,mat_id)
     SG_Eta = eos_stat(3,mat_id)
     SG_Cv = eos_stat(4,mat_id)
     
     rho_min = eos_stat(17,mat_id)
     p_min = eos_stat(18,mat_id)
     ine_min = eos_stat(19,mat_id)

     if(dd.lt.rho_min .or. pp.lt.p_min) then
        print *, " in  eos_ienergy: "
        print *, " rho_min, dd = ", rho_min, dd
        print *, " p_min, pp = ", p_min, pp
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if
     
     ee =  (pp + (SG_Gam-1.d0)*SG_Eta*dd + SG_Gam*SG_Pi)/(SG_Gam-1.d0)/dd

     if(ee.lt.ine_min ) then
        print *, " in  eos_ienergy: "
        print *, " ine_min = ", ine_min
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if

  case(2)

     MG_Gam = eos_stat(1,mat_id)
     MG_rho = eos_stat(2,mat_id)
     MG_c = eos_stat(3,mat_id)
     MG_s = eos_stat(4,mat_id)
     MG_R = eos_stat(5,mat_id)

     rho_min = eos_stat(17,mat_id)
     p_min = eos_stat(18,mat_id)
     ine_min = eos_stat(19,mat_id)
     c_min = eos_stat(20,mat_id)     

     if(dd.lt.rho_min .or. pp.lt.p_min) then
        print *, " in  eos_press_sound: "
        print *, " rho_min, dd = ", rho_min, dd
        print *, " p_min, pp = ", p_min, pp
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if

     x_loc= 1.d0- MG_Rho/dd
     if(1.d0-MG_s*x_loc.le.0.3d0) then
         print *, " + in eos_ienergy: 1.d0-MG_s*x_loc is small = ", 1.d0-MG_s*x_loc
         pause
     endif

     if(x_loc.ge.0) then
        ff= MG_Rho*MG_c**2.d0*x_loc*(1.d0-MG_Gam*x_loc/2.d0)/(1.d0-MG_s*x_loc)**2.d0
     else
        ff= MG_c**2.d0*(dd- MG_Rho) 
     endif
     
     ee = (pp-ff)/MG_Gam/MG_Rho

     if(ee.lt.ine_min ) then
        print *, " in  eos_ienergy: "
        print *, " ine_min = ", ine_min
        print *, " mat_id, eos_type = ", mat_id, eos_type
        continue
        pause
     end if
  
  case default
        
     print *, " in eos_ienergy: this is case of eos_type is not included ."
     pause
     stop

  end select


end subroutine eos_ienergy


