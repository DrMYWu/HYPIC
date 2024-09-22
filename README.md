# HYPIC v1.1
Here is the revision log.
Compared to HYPIC v1.0, this release (v1.1) has the following upgrades:
(1)The electron temperature evolution equation is introduced to describe the electron temperature more reasonably. 
The spatial transport of electron temperature and the heating effect of radio frequency waves on electrons are newly considered.
ref[W. Mingyang, X. Chijie, W. Xiaogang, et al, Plasma Science and Technology 24, 055002 (2022).]
(2)Added the power-mode antenna. User can set the total power of the antenna equal to the sum of the power 
absorbed by the plasma (measured as the power lost by the particles from the boundary) and the Joule heat.
(3)Multiple types of ions can be used at the same case.
(4)0-2 RF frequencies can be used at the same case.
(5)The axial origin (z=0) was moved from the left end of the magnetic mirror to the center of the magnetic mirror.
(6)Optimized the antenna current mode.
(7)Automatic selection of the appropriate time step, according to the nature of energy conservation (loss <1%) in the RK4 format.
(8)Optimized the plotting code (in matlab language).
These improvements increase the computational effort, which can be balanced by appropriately reducing the frequency and particle count.

Added by Mingyang Wu, ymwu@pku.edu.cn or 5927754972@qq.com, Peking University
2024-09-05 20:15


# HYPIC v1.0
HYPIC: A fast hybrid EM PIC-MCC code for ion cyclotron resonance energization in cylindrical coordinate system

Simulation model: 
(1)Maxwell's equations in the frequency domain
(2)momentum equation of ions
(3)The adiabatic electron approximation
(4)Monte-Carlo collision (MCC)

(r,phi,z), axisymmetry df/dphi=0

Support antennas: single loop antenna,single-loop antennas, helical antennas, half-turn antennas, and Nagoya antennas.
 All antennas can be used to study antenna-plasma interactions, but only single loop antenna is available to ICRE simulation.
Support frequency: radio frequency (100k ~ 500M)Hz
Support configurations: linear devices, such as high-power electric propulsion, magnetic mirror, and field-reversed-configuration (FRC), etc.
Support background magnetic field: uniform or generated by coils

Cite: M. Wu, A. Xu, C. Xiao, HYPIC: A fast hybrid EM PIC-MCC code for ion cyclotron resonance energization in cylindrical coordinate system, 
Computer Physics Communications 301, 109207 (2024).

At present, the HYPIC code has been tested in Window 10 and Linux under Visiual Studio 2019 and ifort compiler, with the Intel Math Kernel Library (Parallel) included.

Files structures:

- hypic.f90    ! main program, containing global variables and the main structure of the entire program

./modules         ! sub files of the program
  - fun_ini.f90   !Initialize all parameters such as plasma density, electron temperature, time step, total run time, input and output parameters, etc. 
  - fun_grid.f90 ! Define configuration, simulation domain, array dimensions, etc.
  - fun_b0.f90   ! Set the coil parameters and calculate the background magnetic field
  - fun_irf_antenna.f90 ! Set the antenna type and position
  - fun_fdfd.f90 ! Solve the Maxwell's equations in the frequency domain
  - fun_particles.f90 ! Particle initialization, injection, motion, escape, interpolation, etc.
  - fun_mcc.f90  ! Solve the MCC, including ion-ion (i-i), ion-electron (i-e), electron-electron (e-e), and electron-ion (e-i) collisions.
  - fun_record_display.f90 ! Output running status to the screen and record data to files.


How to run the code:
1. Set configuration and parameters in subroutines and global variables

2. Compile.
   Visiual Studio in Windows 10:
     (1) include all the f90 files to 'Source Files'
     (2) set 'Project'->'Properties'->'Fortran'->'Libraries'->'Use Intel Math Kernel Library'->'Parallel(/Qmkl:parallel)'
     (3) run 'start'
   ifort compiler in linux
     (1) compile with order: 
	     ifort hypic.f90 fun_b0.f90 fun_record_display.f90 fun_grid.f90 fun_ini.f90 fun_irf_antenna.f90 fun_particles.f90 fun_mcc.f90 fun_fdfd.f90 -mkl -o 3-1.x
     (2) run order
	     #insure iswitch_display=  1, display the running state on screen.
	     ./3-1.x
		 or
		 #insure iswitch_display=  0, display nothing when running, suitable for backstage running.
		 ./3-1.x & 
	 
3. Plot. Copy the outputted data into the "data_plot" folder and use the MATLAB program (run the " .m" files) to draw the graphs.
 There is no need to wait for the program to finish running to draw diagrams.


If you meet any problems, please contact:
Mingyang Wu, ymwu@pku.edu.cn or 5927754972@qq.com, Peking University

2024-01-04 13:55


