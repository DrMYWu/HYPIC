    subroutine mcc_main
    use the_whole_varibles
    implicit none
    integer i_c_type

    call find_func_cputime_1_of_2  !-------------------1/2
    call mcc_tau


    i_c_type=1 !i-i collision
    if(abs(t-time_to_run_mcc(i_c_type))<1.001*dt .or. t<1.5*dt )then
        time_to_run_mcc(i_c_type)=time_to_run_mcc(i_c_type)+dt_mcc(i_c_type)
        call mcc_ii(dt_mcc(i_c_type))  !i-i collision
    endif


    i_c_type=2 !e-e collision
    if( i_switch_MHD_Te/=1 .and.  (abs(t-time_to_run_mcc(i_c_type))<1.001*dt .or. t<1.5*dt) )then
        time_to_run_mcc(i_c_type)=time_to_run_mcc(i_c_type)+dt_mcc(i_c_type)
        call mcc_ee(dt_mcc(i_c_type)) !e-e collision
    endif

    i_c_type=3 !i-e collision
    if(abs(t-time_to_run_mcc(i_c_type))<1.001*dt .or. t<1.5*dt )then
        time_to_run_mcc(i_c_type)=time_to_run_mcc(i_c_type)+dt_mcc(i_c_type)
        if(i_switch_MHD_Te/=0)call generate_ve_by_vi_and_Te
        call mcc_ie(dt_mcc(i_c_type))  !i-e collision
    endif

    i_c_type=4 !e-i collision
    if(abs(t-time_to_run_mcc(i_c_type))<1.001*dt .or. t<1.5*dt )then
        time_to_run_mcc(i_c_type)=time_to_run_mcc(i_c_type)+dt_mcc(i_c_type)
        if(i_switch_MHD_Te/=0)call generate_ve_by_vi_and_Te
        call mcc_ei(dt_mcc(i_c_type))  !e-i collision
    endif

    !call find_Qie_by_mcc  !not use

    Ek_i_ave=0.
    do ip=1,np_max
        i_s=np_spe(ip)
        Ek_i_ave=Ek_i_ave+mass_q_i_05_spe(i_s)*sum(v(ip,1:3)**2);
    enddo
    Ek_i_ave=Ek_i_ave/real(np_max) !ev
    
    Ek_e_ave=mass_q_e_05*sum(v_e(1:np_max,1:3)**2)/real(np_max); !ev
    ti_ave=Ek_i_ave/1.5
    te_ave=Ek_e_ave/1.5

    call find_func_cputime_2_of_2(func_time(5))  !-----2/2
    END SUBROUTINE mcc_main



    subroutine mcc_tau
    use the_whole_varibles
    implicit none
    integer i_c_type
    real*8 tau_all(4),tau_min
    real*8 qa,qb,ma,mb,g2_ave,b0_ave,lambda_D,t_ave,coeff_ne_qe_lnA

    ne_mcc=density_ave
    !Ek_i_ave=mass_q_i_05_spe*sum(v(1:np_max,1:3)**2)/real(np_max); !ev
    Ek_i_ave=0.
    do ip=1,np_max
        i_s=np_spe(ip)
        Ek_i_ave=Ek_i_ave+mass_q_i_05_spe(i_s)*sum(v(ip,1:3)**2);
    enddo
    Ek_i_ave=Ek_i_ave/real(np_max) !ev
    
    Ek_e_ave=mass_q_e_05*sum(v_e(1:np_max,1:3)**2)/real(np_max); !ev
    ti_ave=Ek_i_ave/1.5
    te_ave=Ek_e_ave/1.5
    t_ave=.5*(ti_ave+te_ave)

    b0_ave=qe_abs**2/(2.*pi*epson0*3.*t_ave*qe_abs);
    lambda_D=sqrt(epson0*t_ave*qe_abs/(ne_mcc*qe_abs**2));
    ln_A=log(lambda_D/b0_ave)
    coeff_ne_qe_lnA=(ne_mcc*qe_4*ln_A);

    !!when tau_ii>5e-4 s set tau_ii=5e-4 s by changing ne_mcc; to consider the collision effect--start
    !mu_ab=mp/2.;
    !g2_ave=3.*t_ave*qa/mu_ab;
    !b0_ave=qi_spe*qi_spe/(2*pi*epson0*mu_ab*g2_ave);
    !lambda_D=sqrt(epson0*t_ave*qi_spe/(ne_mcc*qi_spe**2));
    !ln_A=log(lambda_D/b0_ave)
    !tau_ii=8*pi*sqrt(2.*ma)*ep_2*((qe_abs*ti_ave)**1.5)/(ne_mcc*qe_4*ln_A);
    !if(tau_ii>5e-4 .and. t>2e-3)then
    !    tp1=tau_ii/5e-4
    !    ne_mcc=density_ave*tp1
    !endif
    !!when tau_ii>5e-4 s set tau_ii=5e-4 s  by changing ne_mcc; to consider the collision effect--end


    ! i-i collision
    ma=mi_spe(1)
    !mb=mp
    tau_ii=8*pi*sqrt(2.*ma)*ep_2*((qe_abs*ti_ave)**1.5)/coeff_ne_qe_lnA;
    !v_ii=1/tau_ii

    ! e-e collision
    ma=me
    !mb=me
    tau_ee=8*pi*sqrt(2.*ma)*ep_2*((qe_abs*te_ave)**1.5)/coeff_ne_qe_lnA;
    !v_ee=1/tau_ee

    ! i-e collision
    ma=mi_spe(1)
    mb=me
    !tau_ie=8*pi*sqrt(2.*ma)*ep_2*((qe_abs*t_ave)**1.5)/(ne_mcc*qe_4*ln_A);
    tau_ie=(4.*pi*epson0)**2*3*ma*((qe_abs*te_ave)**1.5)/(8*sqrt(2*pi*mb)*coeff_ne_qe_lnA);
    !v_ee=1/tau_ee

    ! e-i collision
    ma=me
    mb=mi_spe(1)
    mu_ab=ma*mb/(ma+mb);
    tau_ei=         4.*pi*ep_2*sqrt(mu_ab)*((qe_abs*t_ave)**1.5)/coeff_ne_qe_lnA;
    !tau_ei=sqrt(pi)*4.*pi*ep_2*sqrt(ma)*((qe_abs*t_ave)**1.5)/coeff_ne_qe_lnA;
    !v_ee=1/tau_ee

    tau_all(1)=tau_ii;
    tau_all(2)=tau_ee;
    tau_all(3)=tau_ie;
    tau_all(4)=tau_ei;

    dt_over_dt_mcc(4)=5e-2 !other:1e-2
    dt_mcc(1:4)=tau_all(1:4)*dt_over_dt_mcc(1:4)
    do i_c_type=1,4
        if(dt_mcc(i_c_type)/dt>50.)dt_mcc(i_c_type)=dt*50. !run mover (dt) 50 times and then at least run mcc 1 time
        if(dt_mcc(i_c_type)<dt)dt_mcc(i_c_type)=dt
    enddo

    !tau_min=minval(tau_all(1:3))
    !tp1=5e-2*tau_min
    !if( (tp1*1e-3)>dt)tp1=1e-3*tau_min  !run mover (dt) 50 times and then at least run mcc 1 time
    !if(   tp1     <dt)tp1=dt
    !dt_mcc=tp1
    !!dt_mcc=dt  !test!!
    continue
    END SUBROUTINE mcc_tau


    subroutine mcc_ii(dt_mcc_value)
    use the_whole_varibles
    implicit none
    integer*4::n_rand(1:np_max),np_d2,i1,i2,itp_3,i_s1,i_s2
    real*8 qa,qb,ma,mb,dt_mcc_value
    real*8 ::vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qi_spe(1);
    qb=qi_spe(1);
    ma=mi_spe(1);
    mb=mi_spe(1);
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc_value*(qa*qb/(epson0*mu_ab))**2;
    if(mod(np_max,2)/=0) then
        np_d2=nint((np_max-1.)/2.)
    else
        np_d2=nint(np_max/2.)
    endif
    call mcc_randperm(np_max,n_rand,0)

    do ip=1,np_d2
        i1=n_rand(ip);
        i_s1=np_spe(i1)
        qa=qi_spe(i_s1)
        ma=mi_spe(i_s1)

        i2=n_rand(ip+np_d2);
        i_s2=np_spe(i2)
        qb=qi_spe(i_s2)
        mb=mi_spe(i_s2)

        vv_inj_bef=v(i1,1:3)
        vv_tar_bef=v(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v(i1,1:3)=vv_inj_aft
        v(i2,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ii


    subroutine mcc_ee(dt_mcc_value)
    use the_whole_varibles
    implicit none
    integer*4::n_rand(1:np_max),np_d2,i1,i2,itp_3
    real*8 qa,qb,ma,mb,dt_mcc_value
    real*8 ::vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qe_abs;
    qb=qe_abs;
    ma=me;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc_value*(qa*qb/(epson0*mu_ab))**2;
    if(mod(np_max,2)/=0) then
        np_d2=nint((np_max-1.)/2.)
    else
        np_d2=nint(np_max/2.)
    endif
    call mcc_randperm(np_max,n_rand,0)

    do ip=1,np_d2
        i1=n_rand(ip);
        i2=n_rand(ip+np_d2);
        vv_inj_bef=v_e(i1,1:3)
        vv_tar_bef=v_e(i2,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v_e(i1,1:3)=vv_inj_aft
        v_e(i2,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ee



    subroutine mcc_ie(dt_mcc_value)
    use the_whole_varibles
    implicit none
    integer*4 n_rand1(1:np_max),n_rand2(1:np_max),np_d2,itp_3
    real*8 qa,qb,ma,mb,dt_mcc_value
    real*8    vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qi_spe(1);
    qb=qe_abs;
    ma=mp;
    mb=me;
    mu_ab=ma*mb/(ma+mb);
    sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc_value*(qa*qb/(epson0*mu_ab))**2;
    !call mcc_randperm(np_max,n_rand1,1)
    !call mcc_randperm(np_max,n_rand2,1)

    do ip=1,np_max
        i_s=np_spe(ip)
        qa=qi_spe(i_s)
        ma=mi_spe(i_s)

        vv_inj_bef=v(ip,1:3)
        vv_tar_bef=v_e(ip,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v(ip,1:3)=vv_inj_aft
        v_e(ip,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ie



    subroutine mcc_ei(dt_mcc_value)
    use the_whole_varibles
    implicit none
    integer*4 n_rand1(1:np_max),n_rand2(1:np_max),np_d2,itp_3
    real*8 qa,qb,ma,mb,v_ref,t_ref,dt_mcc_value
    real*8    vv_inj_bef(3),vv_tar_bef(3),vv_inj_aft(3),vv_tar_aft(3)

    qa=qe_abs;
    qb=qi_spe(1);
    ma=me;
    mb=mp;
    mu_ab=ma*mb/(ma+mb);
    v_ref=sqrt(qa*.5*(ti_ave+te_ave)/mu_ab)
    t_ref=tau_ei;
    sab_part=v_ref**3*dt_mcc_value/t_ref
    !!sab_part=ne_mcc*ln_A/(4*pi)*dt_mcc_value*(qa*qb/(epson0*mu_ab))**2;
    !call mcc_randperm(np_max,n_rand1,1)
    !call mcc_randperm(np_max,n_rand2,1)

    do ip=1,np_max
        i_s=np_spe(ip)
        qb=qi_spe(i_s)
        mb=mi_spe(i_s)

        vv_inj_bef=v_e(ip,1:3)
        vv_tar_bef=v(ip,1:3)
        call mcc_sub(ma,mb,qa,qb,vv_inj_bef,vv_tar_bef,vv_inj_aft,vv_tar_aft)
        v_e(ip,1:3)=vv_inj_aft
        v(ip,1:3)=vv_tar_aft
    enddo
    endsubroutine mcc_ei


    subroutine mcc_randperm(ncc,n_rand,itype)
    use the_whole_varibles
    implicit none
    integer*4::ncc,ncc_d2,ilen,n_len,itp2,itype
    integer*4::n_bef(1:ncc),n_rand(1:ncc)
    real*8::rand_0_1
    integer*4 itp0

    if(itype==0) then
        ! i-i or e-e grouping
        if(mod(ncc,2)/=0) then
            n_len=ncc-1
        else
            n_len=ncc
        endif
        ncc_d2=n_len/2

        do itp0=1,n_len
            n_bef(itp0)=itp0
        enddo
        n_rand(1:n_len)=0
        do itp0=1,n_len
            ilen=n_len-itp0+1;
            call random_number(rand_0_1)
            itp2=floor(ilen*rand_0_1)+1;
            n_rand(itp0)=n_bef(itp2);
            n_bef(itp2)=n_bef(ilen);
        enddo
        n_rand(ncc)=ncc

    else

        ! i-e or e-i grouping
        n_len=ncc;
        do itp0=1,n_len
            n_bef(itp0)=itp0
        enddo
        n_rand(1:n_len)=0
        do itp0=1,n_len
            ilen=n_len-itp0+1;
            call random_number(rand_0_1)
            itp2=floor(ilen*rand_0_1)+1;
            n_rand(itp0)=n_bef(itp2);
            n_bef(itp2)=n_bef(ilen);
        enddo

    endif
    END SUBROUTINE



    subroutine mcc_sub(ma,mb,qa,qb,v_inj_bef,v_tar_bef,v_inj_aft,v_tar_aft)
    use the_whole_varibles
    implicit none
    Real*8 ::ma,mb,qa,qb
    real*8 ::v_inj_bef(3),v_tar_bef(3),v_inj_aft(3),v_tar_aft(3)
    !Real*8 ::ne_mcc,mu_ab
    Real*8 ::rand_0_1,s_ab,a_ab
    Real*8 ::g_ab,ep_ab,cos_ep,sin_ep,g_perp,cosx,sinx
    Real*8 ::tp_min,u_ab
    Real*8 ::tp_ab(1:3)
    real*8 ::g3_ab(1:3),h3_ab(1:3)
    integer*4 itp0

    g3_ab(1:3)=v_inj_bef(1:3)-v_tar_bef(1:3)+1e-5; !+1e-5 to avoid a3_ab=0  , very important!
    !need know g_ab=va-vb
    g_ab=sqrt(g3_ab(1)**2+g3_ab(2)**2+g3_ab(3)**2);
    s_ab=sab_part*g_ab**(-3); !find s

    !------find A--------------!
    if (s_ab<0.01) then
        a_ab=1/s_ab;
    elseif (s_ab<4) then
        !            [tp_min,itp0]
        tp_min=minval(abs(s_ab-s_tab));
        itp0=minloc(abs(s_ab-s_tab),1)
        if (s_ab-s_tab(itp0)>0)then
            tp1=s_ab-s_tab(itp0);
            tp2=s_tab(itp0+1)-s_ab;
            a_ab=(tp1*a_tab(itp0+1)+tp2*a_tab(itp0))/(tp1+tp2);
        else
            tp1=s_tab(itp0)-s_ab;
            tp2=s_ab-s_tab(itp0-1);
            a_ab=(tp1*a_tab(itp0-1)+tp2*a_tab(itp0))/(tp1+tp2);
        endif
    else
        a_ab=3*exp(-s_ab);
    endif
    !------find A--------------!

    !find cosx
    call random_number(rand_0_1)
    u_ab=rand_0_1
    if (s_ab<0.01) then
        cosx=1+s_ab*log(u_ab);
    elseif (s_ab<6) then
        cosx=1/a_ab*log(exp(-a_ab)+2*u_ab*sinh(a_ab));
    else
        cosx=2*u_ab-1;
    endif
    sinx=sqrt(1-cosx**2);

    !find v after collision
    call random_number(rand_0_1)
    ep_ab=2*pi*rand_0_1
    cos_ep=cos(ep_ab);
    sin_ep=sin(ep_ab);
    g_perp=sqrt(g3_ab(2)**2+g3_ab(3)**2);
    h3_ab(1)=g_perp*cos_ep;
    h3_ab(2)=-(g3_ab(2)*g3_ab(1)*cos_ep+g_ab*g3_ab(3)*sin_ep)/g_perp;
    h3_ab(3)=-(g3_ab(3)*g3_ab(1)*cos_ep-g_ab*g3_ab(2)*sin_ep)/g_perp;
    tp_ab(1:3)=(g3_ab*(1-cosx)+h3_ab*sinx);

    v_inj_aft=v_inj_bef-mb/(ma+mb)*tp_ab;
    v_tar_aft=v_tar_bef+ma/(ma+mb)*tp_ab;
    END SUBROUTINE mcc_sub



    subroutine mcc_constant
    use the_whole_varibles
    implicit none
    s_tab(1:22)=(/0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,  &
        &    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 2., 3., 4./);
    a_tab(1:22)=(/100.5, 50.5, 33.84, 25.50, 20.50, 17.17, 14.79, 13.01, 11.62, 10.51, 5.516,  &
        &    3.845, 2.987, 2.448, 2.067, 1.779, 1.551, 1.363, 1.207, 0.4105, 0.1496 ,0.05496/);
    END SUBROUTINE



    SUBROUTINE generate_ve_by_vi_and_Te
    USE the_whole_varibles
    implicit none
    integer::I_r,I_z,coeff_vth
    real*8::rand_0_1,rtp
    real*8::r_p(1:np_max)

    coeff_vth=sqrt(qe_abs/me)

    r_p(1:np_max)=sqrt(x(1:np_max,1)**2+x(1:np_max,2)**2)
    ir1_iz1_grid(1:np_max,1)=int((r_p(1:np_max)-r(1))/dr)+1
    ir1_iz1_grid(1:np_max,2)=int((x(1:np_max,3)-z(1))/dz)+1

    do ip=1,np_max
        ir=ir1_iz1_grid(ip,1)
        iz=ir1_iz1_grid(ip,2)
        ve_ex=coeff_vth*sqrt(te_mhd(ir,iz))

        call random_number(tp1)
        call random_number(tp2)
        tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
        v_e(ip,1)=ve_ex*tp3

        call random_number(tp1)
        call random_number(tp2)
        tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
        v_e(ip,2)=ve_ex*tp3

        call random_number(tp1)
        call random_number(tp2)
        tp3=sqrt(-2*log(tp1))*cos(2*pi*tp2)
        v_e(ip,3)=ve_ex*tp3
    end do
    !Ek_i(1:np_max)=sqrt(v(1:np_max,1)**2+v(1:np_max,2)**2+v(1:np_max,3)**2)
    continue
    END SUBROUTINE generate_ve_by_vi_and_Te

    !
    !SUBROUTINE find_Qie_by_mcc
    !USE the_whole_varibles
    !implicit none
    !Q_ie=0
    !
    !do ip=1,np_max
    !    ir=ir1_iz1_grid(ip,1)
    !    iz=ir1_iz1_grid(ip,2)
    !    Q_ie(ir,iz)=Q_ie(ir,iz)+Ek_i(ip)-(sqrt(v(ip,1)**2+v(ip,2)**2+v(ip,3)**2))
    !end do
    !Q_ie=mass_q_i_05_spe*Q_ie/minval(dt_mcc)
    !END SUBROUTINE find_Qie_by_mcc