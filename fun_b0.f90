    subroutine background_b0
    use the_whole_varibles
    implicit none
    integer::mr,mz,N_coils,i_coil,N_fermag,i_fermag !Ferromagnetic media
    integer::ir_tp1,iz_tp1,ir_tp2,iz_tp2,ite
    integer::Nr1_b0,Nr2_b0,Nz1_b0,Nz2_b0,iswitch_fermag
    integer*4::ite_max
    integer::I_r,I_z,ir1,ir2,iz1,iz2
    real*8   ::rd_b0,zd_b0,pre,db,h_b0,dr_b0,dz_b0
    real*8   ::a1,a2,a3,a4,a5,a6,ztp,rtp,s_drz_b0,s1,s2,s3,s4,min_r,min_z,dr01,dr02,dz01,dz02
    real*8::dr2_b0,dz2_b0,rs_b0,zs_b0,d2r_b0,d2z_b0,z_magcenter
    real*8::tp1_b0,tp2_b0,tp3_b0,tp4_b0,mu_tp,br_tp,bz_tp,b0_tp,tp5_b0!,tp6_b0,tp7_b0,tp8_b0,tp9_b0
    real*8, allocatable :: I_th(:),j_th(:),posi_coils(:,:),posi_fermag(:,:),Ath_ite(:)
    integer, allocatable :: Nrz_posi(:,:),Nrz_fermag(:,:),sign_fermag(:,:)
    real*8, allocatable :: B_c(:,:,:),jc(:,:),A_th(:,:),A_old(:,:),mu_rel(:,:),jm_th(:,:)
    real*8, allocatable :: r_b0(:),r2_b0(:),z_b0(:)
    real*8  leng_z_coil_1,leng_z_coil_2,leng_z_coil_3

    !background magnetic field generated by coils
    !here (subroutine background_b0) have an independent (denser) grid
    !different from the sparse grid in main function
    !i.e.:
    !denser grid : fun_b0
    !sparser grid: main function

    rs_b0=0.    !-calulation region
    rd_b0=0.5;
    zs_b0=-1.
    zd_b0=1.; !-calulation region
    !z_magcenter=zs_b0+(zd_b0-zs_b0)*0.5

    N_coils=3;  !number of the coil
    N_fermag=1; !number of the ferromagnetic medium
    pre=1e-14 ; !Iteration Accuracy
    ite_max=10000000

    mz=601
    !if(iswtich_inner_dipole==1)mz=401

    dz_b0=(zd_b0-zs_b0)/(real(mz)-1.);
    h_b0=dz_b0
    dr_b0=dz_b0
    mr=nint((rd_b0-rs_b0)/dr_b0+1);
    !mr=nint(nr*real((rd_b0-rs_b0)/rl));
    if(iswitch_display==1)then
        !write(*,*)0.65
        write(*,*)'the background_B0 subroutine is running ......'
        write(*,*)'mr,mz=',mr,mz
    endif

    allocate (I_th(0:N_coils),j_th(0:N_coils))
    allocate (posi_coils(1:4,0:N_coils),Nrz_posi(1:4,0:N_coils))
    allocate (B_c(1:mr,1:mz,1:2),jc(1:mr,1:mz),A_th(1:mr,1:mz),A_old(1:mr,1:mz)) !B_c contains only Br and Bz
    allocate (r_b0(1:mr),r2_b0(1:mr),z_b0(1:mz),mu_rel(1:mr,1:mz),sign_fermag(1:mr,1:mz),jm_th(1:mr,1:mz))
    allocate (posi_fermag(1:4,1:N_fermag),Nrz_fermag(1:4,1:N_fermag),Ath_ite(1:ite_max))
    jm_th=0.;    B_c=0;    jc=0;
    A_th=0;    A_old=0; Ath_ite=0;
    sign_fermag=0;    mu_rel=1.;
    posi_coils=0.
    I_th=0.
    
    
    !The ampere-turns and locations of the individual coils are given below.
    leng_z_coil_1=0.6   !m
    leng_z_coil_2=0.05  !m
    leng_z_coil_3=0.05  !m
 
    i_coil=1; I_th(i_coil)=1.0*leng_z_coil_1*1e5;!A 
    i_coil=2; I_th(i_coil)=3.5*leng_z_coil_2*1e5;!A
    i_coil=3; I_th(i_coil)=3.5*leng_z_coil_3*1e5;!A 
    
    !(1:4)=(r1,r2,z1,z2)
    ! coil 1
    i_coil=1;
    posi_coils(1,i_coil)=0.1;
    posi_coils(2,i_coil)=posi_coils(1,i_coil)+0.02;
    posi_coils(3,i_coil)=-0.3;
    posi_coils(4,i_coil)=posi_coils(3,i_coil)+leng_z_coil_1;
    
    ! coil 2
    i_coil=2;
    posi_coils(1,i_coil)=0.07;
    posi_coils(2,i_coil)=0.1;
    posi_coils(3,i_coil)=-0.35;
    posi_coils(4,i_coil)=posi_coils(3,i_coil)+leng_z_coil_2;

    ! coil 3
    i_coil=3;
    posi_coils(1,i_coil)=0.07;
    posi_coils(2,i_coil)=0.1;
    posi_coils(3,i_coil)=0.3;
    posi_coils(4,i_coil)=posi_coils(3,i_coil)+leng_z_coil_3;
    
    
    !(1:4)=(r1,r2,z1,z2)
    ! Ferromagnetic media
    iswitch_fermag=0; !1-on, 0-off
    i_fermag=1;
    posi_fermag(1,i_fermag)=0.1;
    posi_fermag(2,i_fermag)=0.15;
    posi_fermag(3,i_fermag)=0.3;!-0.2
    posi_fermag(4,i_fermag)=0.5;!


    do ir=1,mr
        r_b0(ir)=rs_b0+dr_b0*(real(ir)-1.)
    end do
    r2_b0=r_b0-dr_b0*0.5;
    do iz=1,mz
        z_b0(iz)=zs_b0+dz_b0*(real(iz)-1.)
    end do


    dr2_b0=dr_b0*dr_b0;
    dz2_b0=dz_b0*dz_b0;
    d2r_b0=2*dr_b0;
    d2z_b0=2*dz_b0;
    do itp=1,2
        Nrz_posi(itp,:)=nint((posi_coils(itp,:)-rs_b0)/dr_b0);
        Nrz_fermag(itp,:)=nint((posi_fermag(itp,:)-rs_b0)/dr_b0);
    enddo
    do itp=3,4
        Nrz_posi(itp,:)=nint((posi_coils(itp,:)-zs_b0)/dz_b0);
        Nrz_fermag(itp,:)=nint((posi_fermag(itp,:)-zs_b0)/dz_b0);
    enddo

    j_th=I_th/(dr_b0*(Nrz_posi(2,:)-Nrz_posi(1,:)+1.)*dz_b0*(Nrz_posi(4,:)-Nrz_posi(3,:)+1.));

    do i_coil=0,N_coils
        ir_tp1=Nrz_posi(1,i_coil);
        ir_tp2=Nrz_posi(2,i_coil);
        iz_tp1=Nrz_posi(3,i_coil);
        iz_tp2=Nrz_posi(4,i_coil);
        if(iz_tp1<1)iz_tp1=1
        if(iz_tp2<1)iz_tp2=1
        if(ir_tp1<1)ir_tp1=1
        if(ir_tp2<1)ir_tp2=1

        jc(ir_tp1:ir_tp2,iz_tp1:iz_tp2)=j_th(i_coil);
    enddo

    do i_fermag=1,N_fermag
        ir_tp1=Nrz_fermag(1,i_fermag);
        ir_tp2=Nrz_fermag(2,i_fermag);
        iz_tp1=Nrz_fermag(3,i_fermag);
        iz_tp2=Nrz_fermag(4,i_fermag);
        if(iz_tp1<1)iz_tp1=1
        if(iz_tp2<1)iz_tp2=1
        if(ir_tp1<1)ir_tp1=1
        if(ir_tp2<1)ir_tp2=1

        sign_fermag(ir_tp1:ir_tp2,iz_tp1:iz_tp2)=1;
    enddo

    ite=0;
    db=1;
    do while((ite<ite_max) .and. (db>=pre))
        ite=ite+1;
        A_old=A_th;

        do iz =2,mz-1
            do  ir =2,mr-1
                tp1_b0=2/dr2_b0+2/dz2_b0+1/(r_b0(ir)*r_b0(ir));
                tp2_b0=mu0*jc(ir,iz);
                tp3_b0=(A_th(ir+1,iz)+A_th(ir-1,iz))/dr2_b0+(A_th(ir,iz+1)+A_th(ir,iz-1))/dz2_b0;
                tp4_b0=1/r_b0(ir)*(A_th(ir+1,iz)-A_th(ir-1,iz))/d2r_b0;
                A_th(ir,iz)=(tp2_b0+tp3_b0+tp4_b0)/tp1_b0;
            enddo
        enddo

        !ir=1; A_th(ir,2:mz-1)=0 !2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
        !iz=1; A_th(2:mr-1,iz)=A_th(2:mr-1,iz+1);
        !ir=mr; A_th(ir,2:mz-1)=A_th(ir-1,2:mz-1)!0;
        !iz=mz; A_th(2:mr-1,iz)=A_th(2:mr-1,iz-1);

        !ir=1; A_th(ir,2:mz-1)=2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
        !iz=1; A_th(2:mr-1,iz)=A_th(2:mr-1,iz+1);
        !ir=mr; A_th(ir,2:mz-1)=0;
        !iz=mz; A_th(2:mr-1,iz)=A_th(2:mr-1,iz-1);

        ir=1; A_th(ir,2:mz-1)=2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
        iz=1; A_th(2:mr-1,iz)=0;
        ir=mr; A_th(ir,2:mz-1)=0;
        iz=mz; A_th(2:mr-1,iz)=0;

        ir=3;iz=mz/2;
        Ath_ite(ite)=A_th(ir,iz);!maxval()
        db=maxval(abs(A_old-A_th));
    enddo

    ir_tp1=minval(Nrz_fermag(1:2,:));
    ir_tp2=maxval(Nrz_fermag(1:2,:));
    iz_tp1=minval(Nrz_fermag(3:4,:));
    iz_tp2=maxval(Nrz_fermag(3:4,:));

    do ir =ir_tp1,ir_tp2
        do  iz =iz_tp1,iz_tp2
            if(sign_fermag(ir,iz)>0.5)then
                br_tp=-(A_th(ir,iz+1)-A_th(ir,iz-1))/d2z_b0;
                bz_tp=(A_th(ir+1,iz)-A_th(ir-1,iz))/d2r_b0+A_th(ir,iz)/r_b0(ir);
                b0_tp=sqrt(br_tp**2+bz_tp**2)
                call find_relative_mu(b0_tp,mu_tp)
                mu_rel(ir,iz)=mu_tp
            endif
        enddo
    enddo

    if(iswitch_fermag>0.5)then
        ir_tp1=minval(Nrz_fermag(1:2,:));
        ir_tp2=maxval(Nrz_fermag(1:2,:));
        iz_tp1=minval(Nrz_fermag(3:4,:));
        iz_tp2=maxval(Nrz_fermag(3:4,:));

        ite=0;
        db=1;
        !mu_rel=1.
        do while((ite<ite_max) .and. (db>=pre))!ite_max
            ite=ite+1;
            A_old=A_th;

            do iz =2,mz-1
                do  ir =2,mr-1
                    tp1_b0=1./(mu_rel(ir,iz)+mu_rel(ir+1,iz));
                    tp2_b0=(r_b0(ir)/r2_b0(ir))/(mu_rel(ir,iz)+mu_rel(ir,iz+1));
                    tp3_b0=(r_b0(ir)/r2_b0(ir+1))/(mu_rel(ir+1,iz)+mu_rel(ir+1,iz+1));
                    tp4_b0=1./(mu_rel(ir,iz+1)+mu_rel(ir+1,iz+1));

                    a1=tp1_b0+tp2_b0+tp3_b0+tp4_b0; !C_ij
                    a2=(r_b0(ir-1)/r2_b0(ir))/(mu_rel(ir,iz)+mu_rel(ir,iz+1));   !B_ij
                    a3=(r_b0(ir+1)/r2_b0(ir+1))/(mu_rel(ir+1,iz)+mu_rel(ir+1,iz+1)); !D_ij
                    a4=1./(mu_rel(ir,iz)+mu_rel(ir+1,iz)); !A_ij
                    a5=1./(mu_rel(ir,iz+1)+mu_rel(ir+1,iz+1));!E_ij
                    !a6=mu0*h_b0**2/8.*(jc(ir,iz)+jc(ir+1,iz)+jc(ir,iz+1)+jc(ir+1,iz+1));
                    a6=mu0*h_b0**2/2.*jc(ir,iz);

                    A_th(ir,iz)=(a6+a2*A_th(ir-1,iz)+a3*A_th(ir+1,iz)+a4*A_th(ir,iz-1)+a5*A_th(ir,iz+1))/a1
                enddo
            enddo

            !ir=1; A_th(ir,2:mz-1)=0 !2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
            !iz=1; A_th(2:mr-1,iz)=A_th(2:mr-1,iz+1);
            !ir=mr; A_th(ir,2:mz-1)=A_th(ir-1,2:mz-1)!0;
            !iz=mz; A_th(2:mr-1,iz)=A_th(2:mr-1,iz-1);

            !ir=1; A_th(ir,2:mz-1)=2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
            !iz=1; A_th(2:mr-1,iz)=A_th(2:mr-1,iz+1);
            !ir=mr; A_th(ir,2:mz-1)=0;
            !iz=mz; A_th(2:mr-1,iz)=A_th(2:mr-1,iz-1);

            ir=1; A_th(ir,2:mz-1)=2*A_th(ir+1,2:mz-1)-A_th(ir+2,2:mz-1);
            iz=1; A_th(2:mr-1,iz)=0;
            ir=mr; A_th(ir,2:mz-1)=0;
            iz=mz; A_th(2:mr-1,iz)=0;

            ir=3;iz=mz/2;Ath_ite(ite)=maxval(A_th);!(ir,iz)
            db=maxval(abs(A_old-A_th));
        enddo
    endif

    do iz =2,mz-1
        do  ir =2,mr-1
            B_c(ir,iz,1)=-(A_th(ir,iz+1)-A_th(ir,iz-1))/d2z_b0;
            B_c(ir,iz,2)=(A_th(ir+1,iz)-A_th(ir-1,iz))/d2r_b0+A_th(ir,iz)/r_b0(ir);
        enddo
    enddo
    B_c(1,:,:)=B_c(2,:,:);
    B_c(mr,:,:)=B_c(mr-1,:,:);
    B_c(:,mz,:)=B_c(:,mz-1,:);
    B_c(:,1,:)=B_c(:,2,:);

    Nr1_b0=nint((rs-rs_b0)/dr_b0);
    if(Nr1_b0<0.5)Nr1_b0=1
    Nr2_b0=nr-1+Nr1_b0
    Nz1_b0=nint((zs-zs_b0)/dz);
    Nz2_b0=nz-1+Nz1_b0

    b0_DC(:,:,2)=0.0D0

    s_drz_b0=dr_b0*dz_b0;
    !Interpolation from denser to sparser grids
    !denser grid : fun_b0
    !sparser grid: mian function

    !s1-s4, squre of bottom left, bottom right, upper right, upper left
    ! ^
    ! |  Radial coordinates r, to up
    ! -> Axial coordinates  z, to right
    !  |s4| |s3|
    !   x(r,z)
    !  |s1| |s2|
    do ir=1,nr
        do iz=1,nz
            rtp=r(ir);
            ztp=z(iz);
            I_r=minloc(abs(rtp-r_b0),1);
            I_z=minloc(abs(ztp-z_b0),1);
            min_r=rtp-r_b0(I_r);
            min_z=ztp-z_b0(I_z);
            if (min_r<0)then
                ir1=I_r-1;
                ir2=I_r;
            else
                ir1=I_r;
                ir2=I_r+1;
            endif
            if (min_z<0)then
                iz1=I_z-1;
                iz2=I_z;
            else
                iz1=I_z;
                iz2=I_z+1;
            endif

            dr01=rtp-r_b0(ir1);
            dr02=r_b0(ir2)-rtp;
            dz01=ztp-z_b0(iz1);
            dz02=z_b0(iz2)-ztp;
            s1=dr01*dz01/s_drz_b0;
            s2=dr01*dz02/s_drz_b0;
            s3=dr02*dz02/s_drz_b0;
            s4=dr02*dz01/s_drz_b0;

            b0_DC(ir,iz,1)=s3*B_c(ir1,iz1,1)+s1*B_c(ir2,iz2,1)+s2*B_c(ir2,iz1,1)+s4*B_c(ir1,iz2,1);
            b0_DC(ir,iz,3)=s3*B_c(ir1,iz1,2)+s1*B_c(ir2,iz2,2)+s2*B_c(ir2,iz1,2)+s4*B_c(ir1,iz2,2);
        enddo
    enddo

    b0_DC(:,:,4)=sqrt(b0_DC(:,:,1)*b0_DC(:,:,1)+b0_DC(:,:,3)*b0_DC(:,:,3))

30  format(<mr>(e12.4,' '))
31  format(<mz>(e12.4,' '))
    open (unit=101,file='2for_plot_B0.txt',status='unknown',iostat=ierror)
    write (101,30)B_c(:,:,1)
    write (101,30)B_c(:,:,2)
    write (101,30)A_th
    close (101)

    open (unit=103,file='2for_plot_r_B0.txt',status='unknown',iostat=ierror)
    write (103,30)r_b0
    close (103)
    open (unit=104,file='2for_plot_z_B0.txt',status='unknown',iostat=ierror)
    write (104,31)z_b0
    close (104)

    continue

32  format(<nr>(e12.4,' '))
    open (unit=106,file='B0_load.txt',status='unknown',iostat=ierror)
    write (106,32)b0_DC(1:nr,1:nz,1)
    write (106,32)b0_DC(1:nr,1:nz,3)
    close (106)

    open (unit=108,file='2for_plot_Ath_iteration.txt',status='unknown',iostat=ierror)
    write (108,33)Ath_ite(1:ite)
    close (108)

33  format(<ite>(e12.4,' '))
    if(iswitch_display==1)then
        !write(*,*)
        write(*,*)'iteration number is',ite
        write(*,*)'the background_B0 subroutine is finished.'
    endif

    continue
    end subroutine



    subroutine find_relative_mu(b0_tp,mu_tp)
    !use the_whole_varibles
    implicit none
    real*8::b0_tp,mu_tp,dz1,dz2,min_val
    integer::i_minloc,iz1,iz2
    real*8::b_table(1:23),mu_table(1:23)
    data b_table /0,0.045,0.49,0.78,0.99,1.13,1.42,1.57,1.621,1.641,1.656,1.697,1.764,1.801,&
        & 1.838,1.87,2.,2.136,2.185,2.22,2.25,2.355,2.556/
    data mu_table /500.,2249.,12250.,12999.,12375.,11299.,7097.,3925.,2701.,2051.,1655.,848.,&
        & 441.,300.,229.,187.,100.,53.,36.4,27.7,22.5,11.77,6.39/

    if(b0_tp<2.5)then
        i_minloc=minloc(abs(b0_tp-b_table),1);
        min_val=b0_tp-b_table(i_minloc)

        if (min_val<0)then
            iz1=i_minloc-1;
            iz2=i_minloc;
        else
            iz1=i_minloc;
            iz2=i_minloc+1;
        endif

        dz1=(b0_tp-b_table(iz1))/(b_table(iz2)-b_table(iz1));
        dz2=(b_table(iz2)-b0_tp)/(b_table(iz2)-b_table(iz1));

        mu_tp=mu_table(iz1)*dz2+mu_table(iz2)*dz1
    else
        mu_tp=1
    endif
    if(b0_tp<1e-6)mu_tp=1
    continue
    end subroutine




    subroutine load_B0
    use the_whole_varibles
    implicit none
    real*8:: b0_load(1:nr,1:nzx2)
    !load the existing magnetic field data directly.

202 format(<nr>(e12.4,' '))
    open (unit=2002,file='B0_load.txt',status='old',action='read',iostat=ierror)
    read (2002,202)b0_load
    close(2002)

    b0_DC(1:nr,1:nz,1)=b0_load(1:nr,1:nz);
    b0_DC(:,:,2)=0.0D0
    b0_DC(1:nr,1:nz,3)=b0_load(1:nr,nz+1:nzx2);
    b0_DC(:,:,4)=sqrt(b0_DC(:,:,1)**2+b0_DC(:,:,3)**2)

    if(iswitch_display==1)then
        write(*,*)'The background B0  has been loaded.'
    endif
    end subroutine load_B0

