import fload
import os
import numpy as np
import csv

###### python routine to create the input files from the static/PDRs run 
# AG & PL 27/02/2018
# use in /output directory, with stat-g1n* the name of directories containing the outputs of STATIC/PFDR runs
# personal note: g1 means RAD = 1
# input parameters for this routine: Avin
# input files for this routine: mhd_phys.out, h2levels.out, H2_lev.out, species.out, mhd_speci.out
# WARNING: designed for options 'NH2_lev = 150' AND 'H2_out = AD' in the input_mhd.in used to run the STATIC/PDR run
# output files for this routine: h2levels.in, species.in, input_mhd_C.template, input_mhd_J.template

#dirlist = os.popen('ls -1d stat-test').read().split()
#print 'dirlist=', dirlist
dirlist=['../output/example_shock']
Avin= 0.1
RAD = 1

for dir in dirlist:
### first step is to get the physical parameters
# at the position where Av = Avin
# physical parameters are Tn, nH, NH2, NCO
    col=[]
    print 'dir=', dir
    file = dir+'/mhd_phys.out'
    f = open(file, 'r')
    header = f.readline().split()
    f.close()
    iAv=header.index('Av')
    Avs=np.loadtxt(file,skiprows=1,usecols=[iAv])
    irow=abs(Avs-Avin).argmin()
    print 'irow=', irow
    Av=Avs[irow]
    print 'Av=', Av
    inputs=np.loadtxt(file,skiprows=1)
    Tn=inputs[irow,header.index('Tn')]
    nH=inputs[irow,header.index('nH')]
    nH0 = inputs[0,header.index('nH')]
    NH2=inputs[irow,header.index('NH2')]
    NCO=inputs[irow,header.index('NCO')]
    print 'Tn=', Tn
    print 'nH=', nH
    print 'nH0=', nH0
    print 'NH2=', NH2
    print 'NCO=', NCO
### second step is to get the fractional abundance of H2 levels
# at the position where Av = Avin
# and to create the corresponding h2levels.in file 
# we read the names of the levels in h2levels.out
    file = dir+'/h2levels.out'
    h2levs_name=np.loadtxt(file,dtype=str,delimiter='    ',skiprows=11)[:,0]
# we read the values in H2_lev.out (!!! requires to have run the stat in AD mode !!!)
    file = dir+'/H2_lev.out'
    print 'file=', file
    inputs=np.loadtxt(file,skiprows=1)
    print 'inputs=', inputs
    h2levs=inputs[irow,4:]
    h2levs=h2levs/sum(h2levs)
# we report both in h2levels.in, plus a header  
    file = dir+'/h2levels.in'
    f=open(file,'w')
    f.write('!============================='+'\n')
    f.write('! H2* steady state abundances'+'\n')
    f.write('!============================='+'\n')
    f.write('{:<15}{:8.2f}{:<2}'.format('!--- T       = ',Tn,' K')+'\n')
    f.write('{:<15}{:8.2e}{:<5}'.format('!--- nH      = ',nH,' cm-3')+'\n')
    f.write('{:<15}{:8.2f}{:<4}'.format('!--- Av      = ',Av,' mag')+'\n')
    f.write('{:<15}{:8.2f}'.format('!--- RAD     = ',RAD)+'\n')
    f.write('{:<15}{:8.2e}{:<5}'.format('!--- N_H2_0  = ',NH2,' cm-2')+'\n')
    f.write('{:<15}{:8.2e}{:<5}'.format('!--- N_CO_0  = ',NCO,' cm-2')+'\n')
    f.write('{:<15}{:<8}'.format('!--- NH2_lev = ',150)+'\n')
    f.write('!============================='+'\n')  
    for i in range(0,150):
           print 'i=', i
           print 'zob',h2levs_name[i]
           f.write('{:<9}{:<4}{:8.2e}'.format(h2levs_name[i],'    ',h2levs[i])+'\n')
    f.close()
### third step is to get the fractional abundance of species
# at the position where Av = Avin
# and to create the corresponding species.in file   
# we read the global structure of the file in species.out
    file = dir+'/species.out'
# doesn't work   species_struct=np.loadtxt(file,dtype={'names': ('index', 'species_name', 'species_init_fa', 'species_for_enth', 'zob'),'formats': ('i3', 'str26', 'f9.3', 'f10', 'str24')}, skiprows=9)
# doesn't work   species_struct=np.loadtxt(file,dtype={'names': ('index', 'species_name', 'species_init_fa', 'species_for_enth', 'zob'),'formats': ('3i', '26str', '9str', '10f', '24str')}, skiprows=9)
# doesn't work   species_struct=np.loadtxt(file,dtype={'names': ('index', 'species_name', 'species_init_fa', 'species_for_enth', 'zob'),'formats': ('{i3}', '{26}', '{9}', '{10}', '{24}')}, skiprows=9)
# almost works    species_struct=np.loadtxt(file,dtype=[('index', '<i8'), ('species_name', 'S26'), ('species_init_fa', 'S9'), ('species_for_enth', '<f8'), ('zob', 'S24')],skiprows=9)
    species_struct=np.genfromtxt(file,dtype=None,names=['index','species_name','species_init_fa','species_for_enth','zob'], delimiter=[3,26,9,10,24], skip_header=9)
# we read the values of initial fractional abundances in mhd_speci.out
    file = dir+'/mhd_speci.out'
    inputs=np.loadtxt(file,skiprows=2)
    species=inputs[irow,9:-2]
# we report both in species.in, plus a header
    file = dir+'/species.in'
    f=open(file,'w')
    f.write('{:<55}{:8.2f}{:<18}'.format('!---- list of chemical species --- Steady state at T = ',Tn,' K ---------------')+'\n')
    f.write('{:<12}{:8.2e}{:<18}'.format('!---- nH  = ',nH ,' cm-3 ----------')+'\n')
    f.write('{:<12}{:8.2f}{:<18}'.format('!---- Av  = ',Av ,' mag -----------')+'\n')
    f.write('{:<12}{:8.2f}{:<18}'.format('!---- RAD = ',RAD,' ---------------')+'\n')
    f.write('{:<15}{:8.2e}{:<18}'.format('!---- N_H2_0 = ',NH2,' cm-2 ----')+'\n')
    f.write('{:<15}{:8.2e}{:<18}'.format('!---- N_CO_0 = ',NCO,' cm-2 ----')+'\n')
    f.write('!---- WARNING : order = neutrals, species on mantles, ions >0, ions <0 ---------'+'\n')  
    f.write('!---- name, composition, initial density(n/nH), formation enthalpy (kCal/mol) --'+'\n')  
    f.write('!-------------------------------------------------------------------------------'+'\n')  
    print 'len=', len(species)
    for i in range(len(species)):
           print 'i=', i
           f.write('{:>3}{:>26}{:8.3e}{:>10}{:>24}'.format(species_struct[i][0],species_struct[i][1],species[i],species_struct[i][3],species_struct[i][4])+'\n')
    f.close()
### fourth step is to create a input_mhd_C.template
    file = dir+'/input_mhd_C.template'
    f=open(file,'w')
    f.write('!---- input files --------------------------------------------------------------'+'\n')
    f.write('example_shock                            ! modele     : output files radix'+'\n')
    f.write('species.in                               ! specfile   : species file (species list / enthalpies / initial abundances)'+'\n')
    f.write('chemistry.in                             ! chemfile   : chemistry file (reaction list / rates)'+'\n')
    f.write('h2levels.in                              ! h2exfile   : h2* file (if none, population initialized depending on op_H2_in)'+'\n')
    f.write('none                                     ! gridfile   : file containing the grid of position - radiation field'+'\n')
    f.write('!---- shock parameters ---------------------------------------------------------'+'\n')
    f.write('C                                        ! shock_type : \'C\' or \'J\', Isoprotonic or isobaric steady state : \'S1\' or \'S2\', Isoprotonic or isobaric PDR: \'P1\' or \'P2\''+'\n')
    f.write('3                                        ! Nfluids    : 1, 2 ou 3'+'\n')
    f.write('1.00E+00                                 ! Bbeta      : Bfield = Bbeta * sqrt(nH) (micro Gauss)'+'\n')
    f.write('AAAA                                     ! Vs_km      : shock speed (km/s)'+'\n')
    f.write('1.00E+03                                 ! DeltaVmin  : initial Vn - Vi (cm s-1)'+'\n')
    f.write('{:<8.2f}{:<18}'.format(nH0,'                                 ! nH_init    : initial nH = n(H) + 2.0 n(H2) + n(H+) (cm-3)')+'\n')
    f.write('1.00E+01                                 ! Tn         : initial gas temperature (n,i,e) (K)'+'\n')
    f.write('3                                        ! op_H2_in   : initial H2 ortho/para ratio (999.9 -> ETL)'+'\n')
    f.write('!---- environment --------------------------------------------------------------'+'\n')
    f.write('5.00E-17                                 ! Zeta       : cosmic ray ionization rate (s-1)'+'\n')
    f.write('1                                        ! F_ISRF     : radiation field spectrum - 1 = Mathis, 2 = Draine'+'\n')
    f.write('{:<8.2e}{:<18}'.format(RAD,'                                 ! RAD        : radiation field intensity (Habing units)')+'\n')
    f.write('{:<8.2e}{:<18}'.format(Avin,'                                 ! Av0        : initial extinction (magnitudes)')+'\n')
    f.write('1                                        ! F_COUP_RAD : perform a full coupling with radiation field transfer (to compute dissociation rates, desorption, ...) - only if RAD neq 0'+'\n')
    f.write('1                                        ! F_AV       : integrate Av or not'+'\n')
    f.write('1                                        ! F_invAv    : use grain coefficients to compute AV/NH (0) or scale grain coefficient to reproduce inv_Av_fac (1)'+'\n')
    f.write('5.34D-22                                 ! inv_Av_fac : AV/NH (Galaxy : 5.34D-22) - only if F_invAv == 1'+'\n')
    f.write('{:<8.2e}{:<18}'.format(NH2,'                                 ! N_H2_0     : column density of H2 buffer (cm-2)')+'\n')
    f.write('{:<8.2e}{:<18}'.format(NCO,'                                 ! N_CO_0     : column density of CO buffer (cm-2)')+'\n')
    f.write('3.50E+00                                 ! vturb      : turbulent velocity (km s-1, used for Doppler broadening in FGK)'+'\n')
    f.write('!---- grain properties ---------------------------------------------------------'+'\n')
    f.write('1                                        ! F_TGR      : compute grain temperature (1) or keep it constant (0)'+'\n')
    f.write('15                                       ! Tgrains    : initial grain temperature (K)'+'\n')
    f.write('1.00E-06                                 ! amin_mrn   : grain MRN minimum radius (in cm)'+'\n')
    f.write('3.00E-05                                 ! amax_mrn   : grain MRN maximum radius (in cm)'+'\n')
    f.write('3.50E+00                                 ! alph_mrn   : grain MRN index'+'\n')
    f.write('2.00E+00                                 ! rho_grc    : grain core volumic mass (g/cm3)'+'\n')
    f.write('1.00E+00                                 ! rho_grm    : grain mantle volumic mass (g/cm3)'+'\n')
    f.write('!---- excitation & cooling -----------------------------------------------------'+'\n')
    f.write('1                                        ! ieqth      : thermal Balance (1 : solved, 0 : fixed T) - only for \'S\' or \'P\''+'\n')
    f.write('1                                        ! Cool_KN    : Kaufman & Neufeld cooling (1) or analytical formula (0)'+'\n')
    f.write('150                                      ! NH2_lev    : number of H2 levels included'+'\n')
    f.write('200                                      ! NH2_lines_out : maximum number of H2 lines in output file'+'\n')
    f.write('BOTH                                     ! H_H2_flag  : H-H2 collisions : DRF, MM or BOTH !'+'\n')
    f.write('1                                        ! iforH2     : Formation on grain model (1, 2, 3, 4)'+'\n')
    f.write('2                                        ! ikinH2     : Kinetic energy of H2 newly formed (1, 2)'+'\n')
    f.write('1                                        ! pumpH2     : H2 pumping by UV photons taken into account (0, 1)'+'\n')
    f.write('50                                       ! NCO_lev    : number of CO levels included'+'\n')
    f.write('!---- numerical parameters -----------------------------------------------------'+'\n')
    f.write('100000                                    ! Nstep_max  : maximum number of integration steps'+'\n')
    f.write('9.99E+99                                 ! timeJ      : shock age (years)'+'\n')
    f.write('1.00E+08                                 ! duration_max : maximum shock duration (years)'+'\n')
    f.write('1.00E-07                                 ! Eps_V      : precision of computation'+'\n')
    f.write('1.00E+14                                 ! XLL        : caracteristic viscous length (cm)'+'\n')
    f.write('!---- output specifications ----------------------------------------------------'+'\n')
    f.write('1                                        ! F_W_HDF5_STD   : write HDF5  output files'+'\n')
    f.write('0                                        ! F_W_HDF5_CHE   : write HDF5  output files'+'\n')
    f.write('1                                        ! F_W_ASCII  : write ASCII output files'+'\n')
    f.write('10000                                    ! Npthdf5    : maximal number of points in HDF5 files'+'\n')
    f.write('5                                        ! Nstep_w    : number of steps between 2 outputs (for ascii and chemical HDF5 files)'+'\n')
    f.write('FD                                       ! speci_out  : data in mhd_speci.out file - \'AD\' (cm-3), \'CD\' (cm-2) or \'FD\' (n(x)/nH)'+'\n')
    f.write('ln(N/g)                                  ! H2_out     : data in H2_lev.out    file - \'AD\' (cm-3), \'CD\' (cm-2) or \'ln(N/g)\''+'\n')
    f.write('local                                    ! line_out   : data in H2_line.out   file - \'local\' (erg/s/cm3) or \'integrated\' (erg/s/cm2/sr)'+'\n')
    f.write('N                                        ! flag_analysis : Output chemical analysis (dominant reactions) (Y/N)'+'\n')
    f.write('!---- developer options --------------------------------------------------------'+'\n')
    f.write('1                                        ! F_SORT     : sort reactions in increasing order before computing derivatives'+'\n')
    f.write('0                                        ! F_CH       : compute CH velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('0                                        ! F_S        : compute S  velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('0                                        ! F_SH       : compute SH velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('\n')
    f.write('!==============================================================================='+'\n')
    f.write('! additional parameter description'+'\n')
    f.write('!==============================================================================='+'\n')
    f.write('iforH2 = 1                               ! Flag : H2 formation on grains'+'\n')
    f.write('                                         !  -1: formation in the v,J = 0,0 and 0,1 levels only'+'\n')
    f.write('                                         !   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)'+'\n')
    f.write('                                         !   1: Proportional to Boltzman Distrib at 17249 K'+'\n')
    f.write('                                         !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)'+'\n')
    f.write('                                         !   3: v = 6, J = 0,1'+'\n')
    f.write('                                         !   4: fraction = relative populations at t, initialised as H2_lev%density'+'\n')
    f.write('                                         !                 and changed during integration'+'\n')
    f.write('ikinH2 = 2                               ! Flag : H2 formation energy released as kinetic energy'+'\n')
    f.write('                                         !   1: 0.5 * (4.4781 - internal)'+'\n')
    f.write('                                         !   2: Inf(1.4927 eV, 4.4781 - internal)'+'\n')
    f.close()
### fifth step is to create a input_mhd_J.template
    file = dir+'/input_mhd_J.template'
    f=open(file,'w')
    f.write('!---- input files --------------------------------------------------------------'+'\n')
    f.write('example_shock                            ! modele     : output files radix'+'\n')
    f.write('species.in                               ! specfile   : species file (species list / enthalpies / initial abundances)'+'\n')
    f.write('chemistry.in                             ! chemfile   : chemistry file (reaction list / rates)'+'\n')
    f.write('h2levels.in                              ! h2exfile   : h2* file (if none, population initialized depending on op_H2_in)'+'\n')
    f.write('none                                     ! gridfile   : file containing the grid of position - radiation field'+'\n')
    f.write('!---- shock parameters ---------------------------------------------------------'+'\n')
    f.write('J                                        ! shock_type : \'C\' or \'J\', Isoprotonic or isobaric steady state : \'S1\' or \'S2\', Isoprotonic or isobaric PDR: \'P1\' or \'P2\''+'\n')
    f.write('1                                        ! Nfluids    : 1, 2 ou 3'+'\n')
    f.write('1.00E-01                                 ! Bbeta      : Bfield = Bbeta * sqrt(nH) (micro Gauss)'+'\n')
    f.write('AAAA                                     ! Vs_km      : shock speed (km/s)'+'\n')
    f.write('1.00E+03                                 ! DeltaVmin  : initial Vn - Vi (cm s-1)'+'\n')
    f.write('{:<8.2f}{:<18}'.format(nH0,'                                 ! nH_init    : initial nH = n(H) + 2.0 n(H2) + n(H+) (cm-3)')+'\n')
    f.write('1.00E+01                                 ! Tn         : initial gas temperature (n,i,e) (K)'+'\n')
    f.write('3                                        ! op_H2_in   : initial H2 ortho/para ratio (999.9 -> ETL)'+'\n')
    f.write('!---- environment --------------------------------------------------------------'+'\n')
    f.write('5.00E-17                                 ! Zeta       : cosmic ray ionization rate (s-1)'+'\n')
    f.write('1                                        ! F_ISRF     : radiation field spectrum - 1 = Mathis, 2 = Draine'+'\n')
    f.write('{:<8.2e}{:<18}'.format(RAD,'                                 ! RAD        : radiation field intensity (Habing units)')+'\n')
    f.write('{:<8.2e}{:<18}'.format(Avin,'                                 ! Av0        : initial extinction (magnitudes)')+'\n')
    f.write('1                                        ! F_COUP_RAD : perform a full coupling with radiation field transfer (to compute dissociation rates, desorption, ...) - only if RAD neq 0'+'\n')
    f.write('1                                        ! F_AV       : integrate Av or not'+'\n')
    f.write('1                                        ! F_invAv    : use grain coefficients to compute AV/NH (0) or scale grain coefficient to reproduce inv_Av_fac (1)'+'\n')
    f.write('5.34D-22                                 ! inv_Av_fac : AV/NH (Galaxy : 5.34D-22) - only if F_invAv == 1'+'\n')
    f.write('{:<8.2e}{:<18}'.format(NH2,'                                 ! N_H2_0     : column density of H2 buffer (cm-2)')+'\n')
    f.write('{:<8.2e}{:<18}'.format(NCO,'                                 ! N_CO_0     : column density of CO buffer (cm-2)')+'\n')
    f.write('3.50E+00                                 ! vturb      : turbulent velocity (km s-1, used for Doppler broadening in FGK)'+'\n')
    f.write('!---- grain properties ---------------------------------------------------------'+'\n')
    f.write('1                                        ! F_TGR      : compute grain temperature (1) or keep it constant (0)'+'\n')
    f.write('15                                       ! Tgrains    : initial grain temperature (K)'+'\n')
    f.write('1.00E-06                                 ! amin_mrn   : grain MRN minimum radius (in cm)'+'\n')
    f.write('3.00E-05                                 ! amax_mrn   : grain MRN maximum radius (in cm)'+'\n')
    f.write('3.50E+00                                 ! alph_mrn   : grain MRN index'+'\n')
    f.write('2.00E+00                                 ! rho_grc    : grain core volumic mass (g/cm3)'+'\n')
    f.write('1.00E+00                                 ! rho_grm    : grain mantle volumic mass (g/cm3)'+'\n')
    f.write('!---- excitation & cooling -----------------------------------------------------'+'\n')
    f.write('1                                        ! ieqth      : thermal Balance (1 : solved, 0 : fixed T) - only for \'S\' or \'P\''+'\n')
    f.write('1                                        ! Cool_KN    : Kaufman & Neufeld cooling (1) or analytical formula (0)'+'\n')
    f.write('150                                      ! NH2_lev    : number of H2 levels included'+'\n')
    f.write('200                                      ! NH2_lines_out : maximum number of H2 lines in output file'+'\n')
    f.write('BOTH                                     ! H_H2_flag  : H-H2 collisions : DRF, MM or BOTH !'+'\n')
    f.write('1                                        ! iforH2     : Formation on grain model (1, 2, 3, 4)'+'\n')
    f.write('2                                        ! ikinH2     : Kinetic energy of H2 newly formed (1, 2)'+'\n')
    f.write('1                                        ! pumpH2     : H2 pumping by UV photons taken into account (0, 1)'+'\n')
    f.write('50                                       ! NCO_lev    : number of CO levels included'+'\n')
    f.write('!---- numerical parameters -----------------------------------------------------'+'\n')
    f.write('100000                                    ! Nstep_max  : maximum number of integration steps'+'\n')
    f.write('9.99E+99                                 ! timeJ      : shock age (years)'+'\n')
    f.write('1.00E+08                                 ! duration_max : maximum shock duration (years)'+'\n')
    f.write('1.00E-07                                 ! Eps_V      : precision of computation'+'\n')
    f.write('1.00E+14                                 ! XLL        : caracteristic viscous length (cm)'+'\n')
    f.write('!---- output specifications ----------------------------------------------------'+'\n')
    f.write('1                                        ! F_W_HDF5_STD   : write HDF5  output files'+'\n')
    f.write('0                                        ! F_W_HDF5_CHE   : write HDF5  output files'+'\n')
    f.write('1                                        ! F_W_ASCII  : write ASCII output files'+'\n')
    f.write('10000                                    ! Npthdf5    : maximal number of points in HDF5 files'+'\n')
    f.write('5                                        ! Nstep_w    : number of steps between 2 outputs (for ascii and chemical HDF5 files)'+'\n')
    f.write('FD                                       ! speci_out  : data in mhd_speci.out file - \'AD\' (cm-3), \'CD\' (cm-2) or \'FD\' (n(x)/nH)'+'\n')
    f.write('ln(N/g)                                  ! H2_out     : data in H2_lev.out    file - \'AD\' (cm-3), \'CD\' (cm-2) or \'ln(N/g)\''+'\n')
    f.write('local                                    ! line_out   : data in H2_line.out   file - \'local\' (erg/s/cm3) or \'integrated\' (erg/s/cm2/sr)'+'\n')
    f.write('N                                        ! flag_analysis : Output chemical analysis (dominant reactions) (Y/N)'+'\n')
    f.write('!---- developer options --------------------------------------------------------'+'\n')
    f.write('1                                        ! F_SORT     : sort reactions in increasing order before computing derivatives'+'\n')
    f.write('0                                        ! F_CH       : compute CH velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('0                                        ! F_S        : compute S  velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('0                                        ! F_SH       : compute SH velocity (1) or adopt neutral velocity (0)'+'\n')
    f.write('\n')
    f.write('!==============================================================================='+'\n')
    f.write('! additional parameter description'+'\n')
    f.write('!==============================================================================='+'\n')
    f.write('iforH2 = 1                               ! Flag : H2 formation on grains'+'\n')
    f.write('                                         !  -1: formation in the v,J = 0,0 and 0,1 levels only'+'\n')
    f.write('                                         !   0: 1/3 of 4.4781 eV in internal energy (=> 17249 K) (Allen, 1999)'+'\n')
    f.write('                                         !   1: Proportional to Boltzman Distrib at 17249 K'+'\n')
    f.write('                                         !   2: Dissociation limit : v = 14, J = 0,1 (4.4781 eV)'+'\n')
    f.write('                                         !   3: v = 6, J = 0,1'+'\n')
    f.write('                                         !   4: fraction = relative populations at t, initialised as H2_lev%density'+'\n')
    f.write('                                         !                 and changed during integration'+'\n')
    f.write('ikinH2 = 2                               ! Flag : H2 formation energy released as kinetic energy'+'\n')
    f.write('                                         !   1: 0.5 * (4.4781 - internal)'+'\n')
    f.write('                                         !   2: Inf(1.4927 eV, 4.4781 - internal)'+'\n')
    f.close()
#    dummy=raw_input()
