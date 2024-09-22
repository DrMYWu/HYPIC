

    !--------------------subroutine antenna current set-----------------!
    subroutine antenna
    use the_whole_varibles
    implicit none
    real*8:: width_r,width_th,width_z,la,ra_p,z1_m,z2_p,elth,elz
    real*8:: z1,z2,z3,z4,ra1
    real*8:: sq_rz,sq_rth,j_loop,j_helical
    integer::nz1m,nz2p1,nz2p,nz1m1,nrz_dz,nrz_dr,nz1,nz2,nz3,nz4,nra1
    complex*16 ::c1,c2
    real*8 :: j_loop_nor,j_helical_nor !ja_value_nor,
    integer iswitch_right_antenna,nz1_right,nz2_right



    iswitch_right_antenna=0; !1-on  0-off: another antenna (f=frequency) symmetrical to the main antenna.
    !antenna set
    ra=0.08
    za=0.0
    width_r=0.01
    width_z=0.04
    width_th=width_z;
    la=0.4     !---only for helical and Nagoya antenna---!


    if(iswitch_RF2==1 .and. mod(maxwellcount,2)==1)then
        !RF2: f=frequency2
        !ra=0.075
        za=-0.3
        !width_r=0.01
        !width_z=0.03
    endif


    ra1=ra+width_r
    z1=za-la/2;
    z2=z1+la;
    z3=za-width_z/2.;
    z4=za+width_z/2.
    z1_m=z1-width_z
    z2_p=z2+width_z
    !---only for helical antenna---!

    rz_ant_region(1)=ra-width_r
    rz_ant_region(2)=ra+width_r*2
    rz_ant_region(3)=za-width_z;
    rz_ant_region(4)=za+width_z;


    !===================== by Xuad =====================!
    nra=1+nint((nr-1)*(ra-rs)/rl);
    nra1=1+nint((nr-1)*(ra1-rs)/rl)
    !nra1=1+nint((nr-1)*ra_p/rl);

    nza=1+nint((nz-1)*(za-zs)/zl)  ;
    nz1=1+nint((nz-1)*(z1-zs)/zl);
    nz2=1+nint((nz-1)*(z2-zs)/zl);
    nz3=1+nint((nz-1)*(z3-zs)/zl)  ;
    nz4=1+nint((nz-1)*(z4-zs)/zl)
    nz1m=1+nint((nz-1)*(z1_m-zs)/zl);
    nz1m1=nz1-1;
    nz2p=1+nint((nz-1)*(z2_p-zs)/zl);
    nz2p1=nz2+1;
    !nzp1=1+nint((nz-1)*zp1/zl);
    !nzp2=1+nint((nz-1)*zp2/zl)

    nz1_right=nz-nz4
    nz2_right=nz-nz3

    !1:4  r1 r2 z1 z2
    !not use the shield_region
    Nrz_shield_region1(1)=1+nint((nr-1)*((ra-width_r-rs))/rl);
    Nrz_shield_region1(2)=nr;
    Nrz_shield_region1(3)=1+nint((nz-1)*((za-width_z*2-zs))/zl)  ;
    Nrz_shield_region1(4)=1+nint((nz-1)*((za+width_z*3-zs))/zl)  ;

    Nrz_shield_region2(1:2)=Nrz_shield_region1(1:2)
    Nrz_shield_region2(3)=nz-Nrz_shield_region1(4)
    Nrz_shield_region2(4)=nz-Nrz_shield_region1(3)


    !===================== by Xuad =====================!
    elth=pi*ra/sqrt((pi*ra)**2+la**2);
    elz=la/sqrt((pi*ra)**2+la**2);

    nrz_dr=nint(nr*width_r/rl)
    nrz_dz=nint(nz*width_z/zl)

    sq_rz=width_r*width_z
    sq_rth=width_r*width_th

    j_loop_nor=1./sq_rz
    j_helical_nor=1./sq_rth


    im=m*i
    imp=im*pi
    ja_m_nor=0.D0;
    !--------------------------Antenna current set--------------------------!
    !!-----single loop----------!
    if(iswitch_antenna_type==0)then
        if (m==0) ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor;
        if(iswitch_right_antenna==1 .and. m==0) ja_m_nor(nra:nra1,nz1_right:nz2_right,2)=j_loop_nor;
    endif
    !!-----single loop----------!


    !-----half loop antenna----------!
    if(iswitch_antenna_type==2)then
        if (m==0) then
            !----------z=center---------!
            ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor/4
        else
            !----------z=center---------!
            ja_m_nor(nra:nra1,nz3:nz4,2)=j_loop_nor*(1-exp(-imp))/(4*imp) !!!
        endif
    endif
    !-----half loop antenna----------!


    if(iswitch_antenna_type==-1 .or. iswitch_antenna_type==1 .or. iswitch_antenna_type==3)then
        if (m==0) then
            !----------z=z1 &z=z2---------!
            ja_m_nor(nra:nra1,nz1m:nz1m1,2)=0.0D0
            ja_m_nor(nra:nra1,nz2p1:nz2p,2)=0.0D0
            !----------z1<z<z2---------!
            ja_m_nor(nra:nra1,nz1:nz2,2)=0;
            ja_m_nor(nra:nra1,nz1:nz2,3)=0;
        else
            !----------z=z1---------!
            ja_m_nor(nra:nra1,nz1m:nz1m1,2)=j_loop_nor*(1-exp(-imp))/(2*imp)
            !----------z=z2---------!
            if(iswitch_antenna_type /= 3) ja_m_nor(nra:nra1,nz2p1:nz2p,2)=j_loop_nor*(1-exp(-imp))/(2*imp)
            if(iswitch_antenna_type == 3) ja_m_nor(nra:nra1,nz2p1:nz2p,2)=-j_loop_nor*(1-exp(-imp))/(2*imp)

            !----------z1<z<z2---------!
            c1=im*width_th/(2*ra)

            !-----Right helical antenna----------!
            if(iswitch_antenna_type==1)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        c2=j_helical_nor*(1-exp(-imp))*exp(-imp*(z(iz)-z1)/la)*( exp(-c1)-exp(c1)  )/(2*imp)
                        ja_m_nor(ir,iz,2)=elth*c2
                        ja_m_nor(ir,iz,3)=elz*c2
                    enddo
                enddo

                !-----Left helical antenna----------!
            elseif(iswitch_antenna_type==-1)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        c2=j_helical_nor*(1-exp(-imp))*exp(imp*(z(iz)-z1)/la)*( exp(-c1)-exp(c1)  )/(2*imp)
                        ja_m_nor(ir,iz,2)=-elth*c2
                        ja_m_nor(ir,iz,3)=elz*c2
                    enddo
                enddo

                !-----Nagoya Type III antenna----------!
            elseif(iswitch_antenna_type==3)then
                do ir=nra,nra1
                    do iz=nz1,nz2
                        ja_m_nor(ir,iz,3)=j_helical_nor*(1-exp(-imp))*( exp(-c1)-exp(c1)  )/(2*imp)
                    enddo
                enddo
            endif

        endif
    endif

    !ja_m_nor=ja_m_nor*1.0D0
    continue
    !pause
    end subroutine



    subroutine find_irf_now
    use the_whole_varibles
    implicit none
    real*8 p_loss, i_t0,max_d_irf,i_pre,i_min,i_max,t_transport,t_stop_change_Irf,ratio_dIrf

    !constant total power mode
    !! I(t0)^2*R+P_loss(t0)=P_total -> t0 means last time
    !! I(t1)^2*R+P_loss(t1)=P_total -> t1 means this time
    !!P_loss(t1)=P_loss(t0)*I(t1)**2/I(t0)**2
    !!I(t1)=sqrt(P_total/(resistance+P_loss(t0)/I(t0)**2))
    t_transport=2e-4 !s
    t_stop_change_Irf=5e-4 !s
    ratio_dIrf=10.*trf/t_transport
    i_min=0.1 ;!A
    i_max=5000 ;!A
    if(t<t_power_on)then
        i_now=i_min
    else
        i_t0=i_now
        p_loss=sum(power_loss_ave(1:4))
        i_pre=sqrt(power/(resistance+p_loss/i_t0**2))

        !if(p_loss<power )then !chane Irf when xxx.i_pre>i_t0 .and. .and. t<t_stop_change_Irf
        !max_d_irf=ratio_dIrf*i_pre
        max_d_irf=ratio_dIrf*(i_pre-i_t0)
        if(max_d_irf>20.)max_d_irf=20.

        i_now=i_pre
        if((i_now-i_t0)>max_d_irf)i_now=i_t0+max_d_irf
        if(i_now<i_min)i_now=i_min
        if(i_now>i_max)i_now=i_max
        !endif
        
        !When t<td*1/4, adjust the antenna current by power.
        !Then, find irf_mean by averaging irf during td/4*<t<td*1/2
        !Using recursive relationship to find irf_mean: I_n=[(n-1)*I_(n-1)+I_n]/n
        !Last, set a constant irf=irf_mean when t>td*1/2;
        if(t<=td/4.)then
            irf_mean=i_now
        elseif(t>td/4. .and. t<td/2)then            
            n_mean=n_mean+1;
            irf_mean=((n_mean-1)*irf_mean+i_now)/n_mean
        else !t>td*2/3 
            i_now=irf_mean
        endif
    endif

    !constant current mode
    if(i_switch_power_mode==0 ) i_now=irf_set

    !find tau_E by set irf=0
    if(t>t_power_end)then
        i_now=0.
        irf2_set=0.
    endif

    power_Joule=i_now**2*resistance;  ! W
    endsubroutine find_irf_now



    subroutine update_and_calculate
    use the_whole_varibles
    implicit none
    integer i1,i2,i4,ir1,ir2,iz1,iz2
    real*8  Erf_tmp(1:nr,1:nz6),irf_tmp
    complex*16:: Erf_PIC(1:nr,1:nz,1:3)


    if(iswitch_RF2 ==  1  .and. mod(maxwellcount,2)==1)then
        irf_tmp=irf2_set
    else
        irf_tmp=i_now
    endif

    e_AC=irf_tmp*e_AC
    b_AC=irf_tmp*b_AC
    ja_AC=irf_tmp*ja_AC
    jp_AC=irf_tmp*jp_AC
    e_int=irf_tmp*e_int

    !--use this won't repeat running absolutly because of small uncertainties caused by pardiso function in FDFD----!
    !Erf_PIC=e_int

    !--reduce the significance digit, to avoid small uncertainties caused by pardiso function in FDFD --------------!
301 format(<nr>(e12.5,' '))
    open (unit=201,file='Erf_tmp.txt',status='unknown',iostat=ierror)
    do itp=1,3
        write (201,301)real(e_int(1:nr,1:nz,itp))
        write (201,301)imag(e_int(1:nr,1:nz,itp))
    enddo
    close(201)
    open (unit=2001,file='Erf_tmp.txt',status='old',action='read',iostat=ierror)
    read (2001,301)Erf_tmp
    close(2001)
    Erf_PIC=0.
    do itp=1,3
        i1=(itp-1)*(2*nz)+1
        i2=i1+nz-1
        i3=i2+1
        i4=i3+nz-1
        Erf_PIC(1:nr,1:nz,itp)=Erf_tmp(1:nr,i1:i2)+i*Erf_tmp(1:nr,i3:i4)
    enddo

    e_int=Erf_PIC

    if(iswitch_RF2 ==  1  .and. mod(maxwellcount,2)==1)then
        Erf6(:,:,4:6)=Erf_PIC
    else
        Erf6(:,:,1:3)=Erf_PIC
    endif
    
    !!test idea source
    !call test_ideal_Erf

    !ir1=Nrz_shield_region1(1)
    !ir2=Nrz_shield_region1(2)
    !iz1=Nrz_shield_region1(3)
    !iz2=Nrz_shield_region1(4)
    !Erf_PIC(ir1:ir2,iz1:iz2,1:3)=0;
    !
    !ir1=Nrz_shield_region2(1)
    !ir2=Nrz_shield_region2(2)
    !iz1=Nrz_shield_region2(3)
    !iz2=Nrz_shield_region2(4)
    !Erf_PIC(ir1:ir2,iz1:iz2,1:3)=0;

    power_depo_2D=irf_tmp*irf_tmp*power_depo_2D
    power_depo_rthz=irf_tmp*irf_tmp*power_depo_rthz
    ptotm=irf_tmp*irf_tmp*ptotm
    ptotal=irf_tmp*irf_tmp*ptotal
    ptotal_boundary=irf_tmp*irf_tmp*ptotal_boundary
    ptotal_inside=irf_tmp*irf_tmp*ptotal_inside
    ptotal_m(:)=irf_tmp*irf_tmp*ptotal_m(:)

    e_output=irf_tmp*e_output
    Xkz_e=irf_tmp*Xkz_e
    max_eth(:)=irf_tmp*max_eth(:)

    if(t>t_power_on .and. t<t_power_end)then
        state_power_on_off=1.
    else
        state_power_on_off=0.
    endif
    end subroutine update_and_calculate
    
    
    !subroutine test_ideal_Erf
    !use the_whole_varibles
    !implicit none
    !integer i1,i2,i4,ir1,ir2,iz1,iz2
    !real*8  Erf_max,rc1,zc1,zc2,zc3,zc4,rw_Erf,zw_Erf,zw_Erf2
    !
    !Erf_max=2e3/zd*0.6;
    !rc1=0.;
    !zc1=-0.;
    !zc2=zs+zd-zc1;
    !rw_Erf=rd;
    !zw_Erf=zd
    !
    !zc3=zs+0.2;
    !zc4=zd-0.2;    
    !zw_Erf2=0.3
    !
    !Erf6=0;
    !do ir=1,nr
    !    do iz=1,nz
    !
    !        tp1=exp(-(r(ir)-rc1)**2/(2*rw_Erf**2));
    !        tp2=exp(-(z(iz)-zc1)**2/(2*zw_Erf**2));
    !        if(z(iz)<zs+0.15 .or. z(iz)>zd-0.15)tp2=0.
    !        
    !        !tp3=exp(-(z(iz)-zc2)**2/(2*zw_Erf**2));
    !        Erf6(ir,iz,2)=Erf_max*tp1*tp2!*(tp2+tp3)
    !        
    !        
    !        
    !        !if(z(iz)<zc4 .and. z(iz)>zc3)then              
    !        !    tp2=exp(-(z(iz)-zc3)**2/(2*zw_Erf2**2));
    !        !    !tp2=0. !test!!!!!!!!!!!
    !        !    
    !        !    tp3=exp(-(z(iz)-zc4)**2/(2*zw_Erf2**2));
    !        !    !tp3=0.!test!!!!!!!!
    !        !
    !        !    Erf6(ir,iz,5)=5.*Erf_max*(tp2+tp3) !tp1*
    !        !endif
    !        
    !    enddo
    !enddo
    !
    !endsubroutine test_ideal_Erf
    