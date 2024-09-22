    SUBROUTINE find_te
    USE the_whole_varibles
    implicit none
    real*8::dt_te0,x_k1(1:nrp,1:nz),x_k2(1:nrp,1:nz),x_k3(1:nrp,1:nz),x_k4(1:nrp,1:nz)
    real*8:: x_k(1:nrp,1:nz),x0(1:nrp,1:nz),x_k_out(1:nrp,1:nz),pe(1:nrp,1:nz),dt_te,ni_bd
    integer irun
    call find_func_cputime_1_of_2  !-------------------1/2

    ni(1:nrp,1:nz)=n_macro*(density_2D(1:nrp,1:nz)+1e-2*maxval(density_2D))
    call find_dTe(dt_te0)

    nrun=nint(dt/dt_te0)
    dt_te=dt/nrun
    do irun=1,nrun

        !!4th order RK difference in time
        !x0=ni*te_mhd;
        !x_k=x0;                call te_RK4(x_k,x_k_out);  x_k1=x_k_out;
        !x_k=x0+dt_te/2.*x_k1;  call te_RK4(x_k,x_k_out);  x_k2=x_k_out;
        !x_k=x0+dt_te/2.*x_k2;  call te_RK4(x_k,x_k_out);  x_k3=x_k_out;
        !x_k=x0+dt_te*x_k3;     call te_RK4(x_k,x_k_out);  x_k4=x_k_out;
        !pe(1:nrp,1:nz)=x0+dt_te/6.*(x_k1+2.*x_k2+2.*x_k3+x_k4)

        !1st order Euler explicit difference in time
        x0=ni*te_mhd;
        x_k=x0;             call te_RK4(x_k,x_k_out);  x_k1=x_k_out;
        pe(1:nrp,1:nz)=x0+dt_te*x_k1

        te_mhd(1:nrp,1:nz)=pe(1:nrp,1:nz)/ni(1:nrp,1:nz)
        te_mhd(1,:)=te_mhd(2,:)
        do iz=1,nz
            do ir=1,nrp
                if(te_mhd(ir,iz)<1.)te_mhd(ir,iz)=1.
                if(te_mhd(ir,iz)>1e4)te_mhd(ir,iz)=1e4
            end do
        enddo
    enddo


    te_2D=te_mhd
    call find_func_cputime_2_of_2(func_time(8))  !-----2/2
    continue
    endSUBROUTINE find_te


    subroutine te_RK4(x_k,x_k_out)
    use the_whole_varibles
    implicit none
    real*8::x_k(1:nrp,1:nz),x_k_out(1:nrp,1:nz),te_k(1:nrp,1:nz),divqe(1:nrp,1:nz)
    real*8 :: dter,dtez,term1,term2,dner,dnez
    integer::nb_type
    real*8::ne_bd,te_bd
    integer::ib_type,ir1,ir2,iz1,iz2
    real*8::div_val,r_tp,dr_tp,dz_tp
    real*8:: nu_0(1:nrp,1:nz),ur_tp,uz_tp
    real*8::para_in(1:5),flux_qe(0:nrp,0:nz,1:3),flux_je(0:nrp,0:nz,1:3)

    divqe=0.;
    flux_je=0.;
    flux_qe=0.;
    te_k=abs(x_k/ni)

    do iz=1,nz
        do ir=1,nrp
            if(te_k(ir,iz)<1.)te_k(ir,iz)=1.
            if(te_k(ir,iz)>1e4)te_k(ir,iz)=1e4
        end do
    enddo

    !-----------------------------jr-----------------------------------------------!
    do iz=2,nz-1
        do ir=1,nrp-1
            para_in(1)=0.5*(b0_DC(ir,iz,1)+b0_DC(ir+1,iz,1))
            para_in(2)=0.5*(b0_DC(ir,iz,3)+b0_DC(ir+1,iz,3))
            para_in(3)=0.5*(te_k(ir,iz)+te_k(ir+1,iz))
            para_in(4)=0.5*(ni(ir,iz)+ni(ir+1,iz))
            para_in(5)=0.5*(ek_ion_2D(ir,iz)+ek_ion_2D(ir+1,iz))/1.5
            call MHD_coeff(para_in)        !calculate the coefficient
            !density flux and energy flux
            dner =(ni(ir+1,iz)-ni(ir,iz))/dr
            dter =(te_k(ir+1,iz)-te_k(ir,iz))/dr
            tp1 = (ni(ir,iz)+ni(ir+1,iz))/2.  !the corresponding ni on half grids
            tp2=ni(ir,iz-1)
            tp3=ni(ir,iz+1)
            dnez=(tp3-tp2)/d2z

            term1=-dner*Da_perp*coeff_1_br2        !--u by MHD----------!
            term2=-dnez*Da_perp*mue_i_brz
            flux_je(ir,iz,1)= (term1+term2)  !
            u_mhd(ir,iz,1)= (term1+term2)/tp1
            ur_tp=u_mhd(ir,iz,1)                   !--u by MHD----------!
            !ur_tp=0.5*(u_pic(ir,iz,1)+u_pic(ir+1,iz,1)) !--u by PIC----------!

            term1=-2.5*de*coeff_1_mue2_br2/coeff_1_mue2_b2*tp1*dter
            term2=2.5*ur_tp*tp1*(.5*(te_k(ir+1, iz)+te_k( ir,iz)))
            flux_qe(ir,iz,1)= term1+term2
        end do
    end do

    !-----------------------------jz----------------------------------------!
    do iz=1,nz-1
        do ir=2,nrp-1
            para_in(1)=0.5*(b0_DC(ir,iz,1)+b0_DC(ir,iz+1,1))
            para_in(2)=0.5*(b0_DC(ir,iz,3)+b0_DC(ir,iz+1,3))
            para_in(3)=0.5*(te_k(ir,iz)+te_k(ir,iz+1))
            para_in(4)=0.5*(ni(ir,iz)+ni(ir,iz+1))
            para_in(5)=0.5*(ek_ion_2D(ir,iz)+ek_ion_2D(ir,iz+1))/1.5
            call MHD_coeff(para_in)
            dnez=(ni(ir,iz+1)-ni(ir,iz))/dz
            dtez=(te_k(ir,iz+1)-te_k(ir,iz))/dz
            tp1 = (ni(ir,iz)+ni(ir,iz+1))/2.  !the corresponding ni on half grids
            tp2=ni(ir-1,iz)
            tp3=ni(ir+1,iz)
            dner=(tp3-tp2)/d2r
            !
            term1=-dner*Da_perp*mue_i_brz        !--u by MHD----------!
            term2=-dnez*Da_perp*coeff_1_bz2
            flux_je(ir,iz,3)= (term1+term2)
            u_mhd(ir,iz,2)= (term1+term2)/tp1
            uz_tp=u_mhd(ir,iz,2)                 !--u by MHD----------!
            !uz_tp=0.5*(u_pic(ir,iz,2)+u_pic(ir,iz+1,2)) !--u by PIC----------!

            term1=-2.5*de*coeff_1_mue2_bz2/coeff_1_mue2_b2*tp1*dtez
            term2=2.5*uz_tp*tp1*(.5*(te_k(ir,iz+1)+te_k(ir,iz)))
            flux_qe(ir,iz,3)= term1+term2
        end do
    end do

    !----------------jr jz boundary------------------!
    do ir=1,nrp
        do iz=1,nz
            ib_type=i_plasma_region(ir,iz)!-10
            if((ib_type>1 .and. ib_type<=4) .or. ib_type>=6 )call MHD_flux_boundary(ib_type,te_k,flux_qe) !,x_k,flux_je
        enddo
    enddo
    !!!--------------------find divJe , divJi , divqe start------------------------------!!
    do ir=1,nrp
        do iz=1,nz
            ib_type=i_plasma_region(ir,iz) !-10
            if(ib_type>=1)then
                call MHD_div(ib_type,r_tp,dr_tp,dz_tp)! .and. ib_type<=5
                tp1=   (flux_qe(ir,iz,1)-flux_qe(ir-1,iz,1))/dr_tp
                tp2=.5*(flux_qe(ir,iz,1)+flux_qe(ir-1,iz,1))/r_tp
                tp3=   (flux_qe(ir,iz,3)-flux_qe(ir,iz-1,3))/dz_tp
                divqe(ir,iz)=tp1+tp2+tp3
            endif
        enddo
    enddo


    nu_0(1:nrp,1:nz)=ni(1:nrp,1:nz)*coeff_Te_1/(qe_abs*te_k(1:nrp,1:nz))**1.5;
    Q_ie(1:nrp,1:nz)=(nu_0(1:nrp,1:nz)*coeff_Te_3*(1+coeff_Te_2*Ek_ion_2D(1:nrp,1:nz)/te_k(1:nrp,1:nz))**(-1.5))* &
        & (Ek_ion_2D(1:nrp,1:nz)/1.5-te_k(1:nrp,1:nz))
    x_k_out(1:nrp,1:nz)=2/3.*(-divqe(1:nrp,1:nz)+ni(1:nrp,1:nz)*Q_ie(1:nrp,1:nz)+power_depo_2D(1:nrp,1:nz));
    endsubroutine te_RK4



    subroutine MHD_flux_boundary(nb_type,te_k,flux_qe)
    use the_whole_varibles
    implicit none
    integer::irb,izb,nb_type
    real*8::te_k(1:nrp,1:nz),ur_tp,uz_tp,vte
    real*8::flux_qe(0:nrp,0:nz,1:3)
    !-----------nb_type----------------!
    !11,23: left boundary   |<-plasma
    !12,24: right boundary     plasma->|
    !13: up  boundary     -^plasma
    !14: bottom boundary  _vplasma
    !15: r=0 boundary
    !16: r and z shared boundary

    vte=coeff_Te_4*sqrt(te_k(ir,iz))

    ur_tp=vte
    uz_tp=vte

    !-----------nb_type=11----------------!
    if(nb_type==11 .or. nb_type==23 .or. nb_type==27  )then
        u_mhd(ir,iz-1,2)= -uz_tp
        flux_qe(ir,iz-1,3)=   -uz_tp*ni(ir,iz)*te_k(ir,iz)  !t2.0  ->2.5  or 5

        !-----------nb_type=2----------------!
    elseif(nb_type==12 .or. nb_type==24 .or. nb_type==28   )then
        u_mhd(ir,iz,2)= uz_tp
        flux_qe(ir,iz,3)=   uz_tp*ni(ir,iz)*te_k(ir,iz)  ! 2.0  ->2.5  or 5

        !-----------nb_type=3----------------!
    elseif(nb_type==13)then
        u_mhd(ir,iz,1)= ur_tp
        flux_qe(ir,iz,1)=   ur_tp*ni(ir,iz)*te_k(ir,iz)  !2.0  ->2.5  or 5

        !-----------nb_type=4----------------!
    elseif(nb_type==14)then
        u_mhd(ir-1,iz,1)= -ur_tp
        flux_qe(ir-1,iz,1)=   -ur_tp*ni(ir,iz)*te_k(ir,iz)  !2.0  ->2.5  or 5
    endif
    endsubroutine MHD_flux_boundary




    subroutine MHD_div(nb_type,r_tp,dr_tp,dz_tp)
    use the_whole_varibles
    implicit none
    integer::nb_type !irb,izb,
    real*8::r_tp,dr_tp,dz_tp !,i1,i2,i4div_val,
    !-----------nb_type----------------!
    !1: left boundary   |<-plasma
    !2: right boundary     plasma->|
    !3: up  boundary     -^plasma
    !4: bottom boundary  _vplasma
    !5: r=0 boundary
    !6: r and z shared boundary

    dr_tp=dr;
    dz_tp=dz;
    r_tp=r(ir)

    !-----------nb_type=1 and nb_type=2----------------!
    if(nb_type==11 .or. nb_type==12 )then ! .or. nb_type==23 .or. nb_type==24
        dz_tp=dz_2

        !-----------nb_type=5----------------!
    elseif(nb_type==15 )then
        dr_tp=dr_2;
        r_tp=dr*0.25

        !-----------nb_type=3----------------!
    elseif(nb_type==13)then !.or. nb_type==23 .or. nb_type==24
        dr_tp=dr_2;
        r_tp=r(ir-1)+dr*0.75!r2(ir)

    elseif(nb_type==14)then
        dr_tp=dr_2;
        r_tp=r(ir)+dr*0.25!r2(ir)

    elseif(nb_type==27 .or.nb_type==28  )then
        dz_tp=dz_2
        dr_tp=dr_2;
        r_tp=r(ir)+dr*0.25!r2(ir)

    elseif( nb_type==25 .or. nb_type==26)then
        dr_tp=dr_2;
        dz_tp=dz_2
        r_tp=dr*0.25

    elseif( nb_type==23 .or. nb_type==24)then
        dr_tp=dr_2;
        dz_tp=dz_2
        r_tp=r(ir-1)+dr*0.75
    endif
    endsubroutine MHD_div



    subroutine MHD_coeff(para_in)
    use the_whole_varibles
    implicit none
    real*8::para_in(1:5)
    real*8::br00,bz00,te00,ne00,bt2,vti,vi,vie,vin,vei,ven,ve,vte
    real*8 :: coeff_1_b2,mue_mui,mue2,di,mui,mue
    real*8:: te_tp,ni_tp,ti_tp
    br00=para_in(1)
    bz00=para_in(2)
    te_tp=para_in(3)
    ni_tp=para_in(4)
    ti_tp=para_in(5)

    bt2=br00**2+bz00**2+1e-10  !B_total^2
    vte=sqrt(8.*te_tp*qe_abs/(pi*me));
    ven=nn*sig_en*vte;
    vei=(te_tp)**(-1.5)*ln_A0*(2.91e-12)*ni_tp;
    ve=ven+vei+1e8; !to avoid De too low
    mue=qe_abs/(me*ve);
    de=te_tp*mue;

    ti_tp=ek_ion_2D(ir,iz)/1.5 !eV
    if(ti_tp<1.)ti_tp=1.
    vti=sqrt(8.*qe_abs*ti_tp/(pi*mp));
    vie=4.80e-14*mass_number(1)**(-0.5)*ti_tp**(-1.5)*ln_A0*ni_tp;
    vin=nn*sig_in*vti;
    vi=vie+vin;
    mui=qe_abs/(mp*vi); !mui
    di=mui*ti_tp    !Di=qe_abs*Ti/(mp*vi);

    mue2=mue*mue
    coeff_1_mue2_br2=1+mue2*br00**2
    coeff_1_mue2_bz2=1+mue2*bz00**2
    coeff_1_mue2_b2=1+mue2*bt2

    mue_mui=mue*mui
    mue_i_brz=mue_mui*br00*bz00
    coeff_1_br2=1+mue_mui*br00**2
    coeff_1_bz2=1+mue_mui*bz00**2
    coeff_1_b2=1+mue_mui*bt2
    Da=(di*mue+de*mui)/(mui+mue) !Da ll
    Da_perp=Da/coeff_1_b2        !Da perp
    endsubroutine MHD_coeff




    subroutine find_dTe(dt_te0)
    use the_whole_varibles
    implicit none
    real*8::para_in(1:5),dTe_2D(nrp,nz),dt_te0
    integer Nte_max
    Nte_max=300
    dTe_2D=sqrt(-1.)
    do iz=1,nz
        do ir=1,nrp
            if(i_plasma_region(ir,iz)>0.5)then
                para_in(1)=b0_DC(ir,iz,1)
                para_in(2)=b0_DC(ir,iz,3)
                para_in(3)=te_mhd(ir,iz)
                para_in(4)=ni(ir,iz)
                para_in(5)=ek_ion_2D(ir,iz)
                call MHD_coeff(para_in)
                dTe_2D(ir,iz)=3.*dr**2/(5.*De)
            endif
        enddo
    enddo

    dt_te0=0.5*minval(dTe_2D)    
    if(dt_te0>dt)dt_te0=dt
    endsubroutine find_dTe


