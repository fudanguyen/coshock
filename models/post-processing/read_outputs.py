from numpy import *
import matplotlib.pyplot as plt
import importlib
import pl; importlib.reload(pl)

option=input('You want to plot the statics? (y/n): ')
if (option=='yes') or (option=='y') or (option=='Y'):
    plt.figure(1)
    #file='../output/'
    file='../output/Cn4v20_steady/mhd_speci.out'
    print ('file name=',file)
    pos = pl.file_load(file)
    dist=pl.v[pos]['distance']
    xCO = pl.v[pos]['x(CO)']
    xH = pl.v[pos]['x(H)']
    plt.loglog(dist,xCO,label='x(CO)')
    plt.loglog(dist,xH,label='x(H)')
    plt.legend()
    
elif (option!='yes') or (option!='y') or (option!='Y'):
    plt.figure(1)
    file='../output/Cn4v20/mhd_phys.out'
    pos = pl.file_load(file)
    dist=pl.v[pos]['distance']
    Tn = pl.v[pos]['Tn']
    Vn = pl.v[pos]['Vn']
    Vi = pl.v[pos]['Vi']
    ti = pl.v[pos]['timeI']
    tn = pl.v[pos]['timeN']
    plt.loglog(dist,Tn,label='Tn')

    plt.figure(2)
    plt.semilogx(dist,Vn/1e5,label='Vn')
    plt.semilogx(dist,Vi/1e5,label='Vi')
    plt.legend()
    
    plt.figure(3)
    plt.loglog(ti,Tn)
    plt.xlabel('$\\rm time_{I}$')
    plt.ylabel('$\\rm Tn $')

    plt.figure(4)
    plt.loglog(ti,ti)
    plt.loglog(ti,tn)
    
else:
    print ('Error!')

plt.show()



