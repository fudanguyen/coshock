
v={}
    
def fl(filename,overwrite=False,nskip=0):
    global v
    from numpy import array
    import re
    import mytools
    if not overwrite:
        v={}
    f=open(filename,'r')
    for i in range(nskip):
        dummy = f.readline()
    lines=f.readlines()
    vars=lines[0].split()
    for var in vars:
        v[var]=[]
    for line in lines[1:]:
        vals=line.split()
        for (var,val) in zip(vars,vals):
            if re.search('[^E]-\d\d\d',val):
                val='0'
            if re.search('NaN',val):
                print 'I found a NaN for var:',var
                val='0'
            if mytools.type_of_value(val)==str:
                val='0'
            v[var].append(float(val))
    for k in v.keys():
        v[k]=array(v[k])


def mp(string):
    from pylab import clf,plot,xlabel,ylabel,loglog,legend,rc
# Thick lines
    rc('lines', linewidth=2)
    (ord,abs)=string.split(';')
    ord=ord.split(',')
    for o in ord:
        plot(v[abs],v[o],label=o)
    legend(loc=0)
    

#fload('../runs/nn2/u5.5.o3/mhd_phys.out')
#fload('../runs/nn2/u5.5.o3/H2_line.out')
#fload('../runs/nn2/u5.5.o3/intensity.out')
#mp('C+(158m),0-0S(1);NH2')
#loglog()

def cs():
    from numpy import reshape
    from pylab import clf,plot,xlabel,ylabel,loglog,legend
    fload('../runs/Cs.dat')
    ch =reshape(v['N(O2)'],(4,6),'F')
    rad=reshape(v['rad'],(4,6),'F')
    av =reshape(v['av'],(4,6),'F')
    clf()
    loglog()
    plot(rad[0][:4],ch[0][:4],'-o',label='Av=0.1')
    plot(rad[0][:4],ch[1][:4],'-o',label='Av=0.3')
    plot(rad[0][:5],ch[2][:5],'-o',label='Av=1')
    plot(rad[0],ch[3],'-o',label='Av=3')
    legend(loc=0)
    xlabel('Radiation field G$_0$')
    ylabel('N(O$^2$) (1/cm$^2$)')
    title("CH$^+$ production in irradiated shocks \n (u=20 km/s, b=1, n$_H$=10$^4$/cm$^3$)")
