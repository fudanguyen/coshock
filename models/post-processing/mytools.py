from vtk import *
import vtk.util.numpy_support
from numpy import prod,size
from numpy import *
 
def cmd(str):
   import os
   print (str)
   os.system(str)

def fact(n):
   f=1
   for i in range(n):
      f=f*(i+1)
   return f

def arr(n,p):
   f=float(n)
   for i in range (p-1):
      f=f*(n-(i+1))
   return f

def farr(n,p):
   f=1.
   for i in range (p-1):
      f=f*((n-(i+1))/n)
   return f

def cnp(n,p):
   return arr(n,p)/fact(p)

def argsh(x):
   import numpy
   return numpy.sign(x)*numpy.log(numpy.sqrt(x*x+1.)+numpy.abs(x))

def argsh10(x):
   import numpy
   return numpy.sign(x)*numpy.log(numpy.sqrt(x*x+1.)+numpy.abs(x))/numpy.log(10.)

# base 10 mantissa and exponent
def frexp_10(decimal):
   import math
   if decimal==0:
      return 0,0
   logdecimal = math.log10(decimal)
   n=int(logdecimal)
   if logdecimal<n:
       n=n-1
   return 10 ** (logdecimal - n), n

# formats a decimal into a string: 3.4167e-24 --> '3.4(-24)'
def myformat(decimal,precision=2):
    # precision is the number of significative digits
    (m,n)=frexp_10(decimal)
    rounded= "%.1f" % m
    return rounded+'('+str(n)+')'
    
def type_of_value(var):
   import ast
   try:
      return type(ast.literal_eval(var))
   except Exception:
      return str

def thick_lines():
   # Thick lines and big fonts
   rc('lines', linewidth=3)
   font = {'family' : 'normal',
           'weight' : 'bold',
           'size'   : 18}
   rc('font', **font)
   rcParams['legend.fontsize']='small'
   
# Read a vtkStrcturedPoints
def read_giorgos_vtk(filename):
   v2n=util.numpy_support.vtk_to_numpy

   Reader = vtkStructuredPointsReader()
   Reader.SetFileName(filename)
   Reader.Update()
   data=Reader.GetOutput()
   dim=data.GetDimensions()
   pd=data.GetPointData()
   dat={}
   for i in range(pd.GetNumberOfArrays()):
      name=pd.GetArrayName(i)
      dat[name]=v2n(data.GetPointData().GetArray(name))
      nspace=size(dat[name])/prod(dim)
      print (name,': nspace=',nspace,' dim=',dim)
      if (nspace>1):
         dat[name]=dat[name].reshape(dim+(nspace,))
      else:
         dat[name]=dat[name].reshape(dim)
   return dat

# Read a vtkRectilinearGrid
def read_dat_vtk(filename):
   v2n=util.numpy_support.vtk_to_numpy

   Reader = vtkRectilinearGridReader()
   Reader.SetFileName(filename)
   Reader.Update()
   data=Reader.GetOutput()
   dim=data.GetDimensions()
   pd=data.GetPointData()
   dat={}
   for i in range(pd.GetNumberOfArrays()):
      name=pd.GetArrayName(i)
      dat[name]=v2n(data.GetPointData().GetArray(name)).reshape(dim)
   return dat

def read_vtk(filename):
   import sys
   import re
   import numpy as np
   f=open(filename,'r')   
   line=f.readline().split()+[0] 
   print (line)
   while(line[0] != 'DATASET'):
      line=f.readline().split()+[0]
   print (line)
   while(line[0] != 'DIMENSIONS'):
      line=f.readline().split()+[0]
   dims=map(int,line[1:4])
   print ('Dimensions:',dims)
   # Read coordinates
   coords=[]
   for ndim in dims:
      line=f.readline().split()
      print (line[0])
      if (int(line[1]) != ndim):
             print (line[1],' is different from ndim=',ndim)
             sys.exit()
      nvalues=0
      x=[]
      while(nvalues<ndim):         
         vals=map(float,f.readline().split())
         nvalues=nvalues+len(vals)
         x=x+vals
      coords.append(x)

   # Read data on the grid
   line=f.readline().split()
   print (line[0])
   if (int(line[1]) != np.prod(dims)):
      print (line[1],' is different from prod(dims)=',np.prod(dims))
      sys.exit()
   line=f.readline().split()
   nfields=int(line[2])
   print (line[0],' nfields=',nfields)
   fields={}
   print ('reading lines..')
   lines=f.readlines()
   print ('..done')
   for field in range(nfields):
      line=lines.pop(0).split()
      print (line)
      name=line[0]
      ndat=int(line[2])
      if (ndat != np.prod(dims)):
         print ('ndat=',ndat,' is different from prod(dims)=',np.prod(dims))
         sys.exit()
      x=np.zeros(ndat)
      nvalues=0
      while(nvalues<ndat):         
         vals=map(float,lines.pop(0).split())
         nval=len(vals)
         x[nvalues:nvalues+nval]=vals
         nvalues=nvalues+len(vals)
      fields[name]=np.reshape(np.array(x),dims)
   f.close()
   return tuple(coords),fields


def read_skel(filename):
   import sys
   import re
   import numpy as np
   f=open(filename+'.a.NDskl','r')   
   line=f.readline().split()   
   while(line[0] != '[FILAMENTS]'):
      line=f.readline().split()
   line=f.readline().split()
   nfil=int(line[0])
   print (nfil,' filaments.')
# Propre:
#   fils=[[]]*nfil
#   for fil in fils:
#      line=f.readline().split()
#      nsamp=int(line[-1])  
#      for j in range(nsamp):
#         coords=map(float,f.readline().split())
#         fil.append(coords)
# Sale:
   fils=[[]]*nfil
   for ifil in range(nfil):
      line=f.readline().split()
      nsamp=int(line[-1])  
      fils[ifil]=[]
      for j in range(nsamp):
         coords=map(float,f.readline().split())
         coords.reverse() # Apparently, coordinates are (y,x)...
         fils[ifil].append(coords)
   return fils

def plot_skel(fils):
   import matplotlib.pyplot as plt
   import numpy as np
   for fil in fils:
      plt.plot(np.array(fil)[:,0],np.array(fil)[:,1],linestyle=':',color='k')
   plt.show()

def merge_skel(fils):
   import sys
   # Get extreme points linked to their filaments
   pts={}
   print (len(fils),' filaments.')
   for ifil in range(len(fils)):
      for extr in [str(fils[ifil][0]), str(fils[ifil][-1])]:
         if (extr in pts.keys()):
            pts[extr].append(ifil)
         else:
            pts[extr]=[ifil]
   print (len(pts),' extreme points.')
   # Merge filaments joined by their points
   # Parse extremeties with exactly two filaments
   todel=[]
   for pt in pts.keys():
       if (len(pts[pt])==2):
           ifila=pts[pt][0]
           ifilb=pts[pt][1]
           if (ifila==ifilb):
   #            print 'Cycle detected'
               break
   #            sys.exit()
           # We want pt in the middle of the merged filament.Merge
           # filament a and b into a. You need to be careful with order..
           if (pt==str(fils[ifilb][-1])):
               fils[ifilb].reverse()
           if (pt==str(fils[ifila][0])):
               fils[ifila].reverse()
           fils[ifila]=fils[ifila]+fils[ifilb]
           # Look at the extremities of filament b:
           pta=str(fils[ifilb][0])
           ptb=str(fils[ifilb][-1])
           # if it is not a cycle:
           if (pta!=ptb):
               # Find the point which is different from pt
               otherpt=ptb
               if (ptb==pt):
                   otherpt=pta
               # For that point, replace filb by fila
               for (i,ifil) in enumerate(pts[otherpt]):
                   if ifil==ifilb:
                       pts[otherpt][i]=ifila
               #   if (ifila in pts[otherpt]):
               #      print 'PROBLEM! I get a cycle...'
               #      sys.exit()
               # delete filb and pt
               # del fils[ifilb]
               todel.append(ifilb)
               del pts[pt] # Remove that point from extremeties
           sys.stdout.write("\r%d  extreme points." % len(pts)) 
   todel.sort()
   todel.reverse()
   for ifil in todel:
       del fils[ifil]
   print ("\n")
   print (len(fils)," filaments.")
   return

def spline_skel(fils,dsmooth):
   import sys
   import scipy.interpolate as interpolate
   import numpy as np

   # define smoothing length
#   xpix=(xx[-1]-xx[0])/len(xx)
#   ypix=(yy[-1]-yy[0])/len(yy)
#   dsmooth=0.5*norm([xpix,ypix])

   # Parse each filament
   for ifil,fil in enumerate(fils):
       if (len(fil)>3):
           # Spline smooth it
           fil2=np.array(fil)+1e-3*dsmooth*np.random.randn(len(fil),2) # Straight lines are no good for fitpack, it seems...
           tck,u=interpolate.splprep(fil2.transpose(),s=len(fil)*(dsmooth)**2)
           sfil=interpolate.splev(np.linspace(0,1.0,len(u)),tck)
           fils[ifil]=np.array(sfil).transpose()

# Find cuts off
def find_cuts(d,fraction):
    from numpy import histogram,array,sqrt,size
    # dissipation histogram
    nbs,bins=histogram(d,bins=sqrt(size(d)))
    # Cumulative histogram of dissipation
    frac=0.0
    s1=[]
    for n,b0,b1 in zip(nbs,bins[:-1],bins[1:]):
        frac=frac+0.5*(b0+b1)*float(n)
        s1.append(frac)
    s1=array(s1)/s1[-1]
    # Center of bins
    xd=array(bins[:-1]+bins[1:])*0.5
    # Extract values between fraction and 1-fraction
    bulkx=xd[(s1>fraction)&(s1<1.-fraction)]
    vmin=bulkx.min()
    vmax=bulkx.max()
    return vmin,vmax

# Build an rgb image from three fields
def rgbim(r,g,b,bounds=[],log=False):
    import Image
    from numpy import sum,shape,reshape,array,zeros,log10,uint8
    etot=r+g+b
    print ('R=',sum(r)/sum(etot),' G=',sum(g)/sum(etot),' B=',sum(b)/sum(etot))
    print ('emin=',etot.min(),' emax=',etot.max(),' bounds=',bounds)
    dim=shape(etot)
    if bounds==[]:
        vmin=etot.min()
        vmax=etot.max()
    else:
        vmin=bounds[0]
        vmax=bounds[1]
    def cuts(v):
        if (v<vmin):
            return 0.
        elif (v>vmax):
            return 1.
        else :
            if log:
                return (log10(v)-log10(vmin))/(log10(vmax)-log10(vmin))
            else:
                return (v-vmin)/(vmax-vmin)
    intensity=array(map(cuts,etot.reshape(dim[0]*dim[1]))).reshape(dim[0:2])
    r=r/etot*intensity
    g=g/etot*intensity
    b=b/etot*intensity
    # cuts off and log scaling for each color.
    rgbArray=zeros(dim[0:2]+(3,))
    maxv=max([r.max(),g.max(),b.max()])
    rgbArray[...,0]=r*255/maxv
    rgbArray[...,1]=g*255/maxv
    rgbArray[...,2]=b*255/maxv
    im=Image.fromarray(uint8(rgbArray))
    return im

def within_pi(a):
   from numpy import pi
#    return a
   return (a+pi)%(2.0*pi)-pi

def within_pi2(a):
   from numpy import pi
#    return a
   return (a+0.5*pi)%(pi)-0.5*pi

# Same for a 2D vector  field
def increments(v_map,lag,norm=2,is_angle=False):
   # Computes the map of the increments of v_map, 
   # abs() for norm!=2, ()**2 for norm=2
   # if is_angle is True, differences are put back within [-pi,pi]
   # also returns the histogram of the increments with their standard deviation.
   # find every relevant shifts and their number
   shifts=[]
   for i in range(-lag-1,lag+2):
      for j in range(-lag,lag+1):
         if ( (lag)**2 <= i*i+j*j < (lag+1)**2 ):
            shifts.append([i,j])
   # for each shift, compute absolute difference(squared ?) between shifted and original array
   # then add it to the CVI map
   nx,ny=shape(v_map)
   iv_map=zeros(shape(v_map))
   ivmax=abs(v_map.max())+abs(v_map.min())
   ivhist=array([0]*nx)
   for shift in shifts:
      inc=v_map-roll(roll(v_map,shift[0],axis=0),shift[1],axis=1)
      if is_angle:
         inc=within_pi2(inc)
      if (norm==2):
         iv_map=iv_map+inc**2
      else:
         iv_map=iv_map+abs(inc)
      if is_angle:
         nbs,bins=histogram(within_pi2(inc),bins=nx,range=[-ivmax,ivmax])
      else:
         nbs,bins=histogram(inc,bins=nx,range=[-ivmax,ivmax])
      ivhist=ivhist+array(nbs)
   cbin=array((bins[:-1]+bins[1:])*0.5)
   sig=sqrt(sum(ivhist*cbin**2)/sum(ivhist))
   iv_map=iv_map/len(shifts)
   if (norm==2):
      iv_map=sqrt(iv_map)
   return iv_map,cbin,double(ivhist),sig


# Same for a 2D vector  field
def increments_2d(u_map,v_map,lag,norm=2,para=True,is_angle=False):
   # Computes the map of the increments of v_map, 
   # abs() for norm!=2, ()**2 for norm=2
   # if is_angle is True, differences are put back within [-pi,pi]
   # also returns the histogram of the increments with their standard deviation.
   # find every relevant shifts and their number
   shifts=[]
   for i in range(-lag-1,lag+2):
      for j in range(-lag,lag+1):
         if ( (lag)**2 <= i*i+j*j < (lag+1)**2 ):
            shifts.append([i,j])
   # for each shift, compute absolute difference(squared ?) between shifted and original array
   # then add it to the CVI map
   nx,ny=shape(v_map)
   iv_map=zeros(shape(v_map[lag:,lag:]))
   ivmax=abs(v_map.max())+abs(v_map.min())
   ivhist=array([0]*nx)
   for shift in shifts:
      incu=u_map-roll(roll(u_map,shift[0],axis=0),shift[1],axis=1)
      incu=incu[lag:,lag:]
      incv=v_map-roll(roll(v_map,shift[0],axis=0),shift[1],axis=1)
      incv=incv[lag:,lag:]
      if para:
         # Parallel increment
         inc=(incu*shift[0]+incv*shift[1])/sqrt(shift[0]**2+shift[1]**2)
      else:
         # Othogonal increment
         inc=(incu*shift[1]+incv*shift[0])/sqrt(shift[0]**2+shift[1]**2)

      if is_angle:
         inc=within_pi2(inc)
      if (norm==2):
         iv_map=iv_map+inc**2
      else:
         iv_map=iv_map+abs(inc)
      if is_angle:
         nbs,bins=histogram(within_pi2(inc),bins=nx,range=[-ivmax,ivmax])
      else:
         nbs,bins=histogram(inc,bins=nx,range=[-ivmax,ivmax])
      ivhist=ivhist+array(nbs)
   cbin=array((bins[:-1]+bins[1:])*0.5)
   sig=sqrt(sum(ivhist*cbin**2)/sum(ivhist))
   iv_map=iv_map/len(shifts)
   if (norm==2):
      iv_map=sqrt(iv_map)
   return iv_map,cbin,double(ivhist),sig


# Same for a 1D field
def increments_1d(v_map,lag,norm=2,is_angle=False):
   # Computes the map of the increments of v_map, 
   # abs() for norm!=2, ()**2 for norm=2
   # if is_angle is True, differences are put back within [-pi,pi]
   # also returns the histogram of the increments with their standard deviation.
   # for each shift, compute absolute difference(squared ?) between shifted and original array
   # then add it to the CVI map
   nx=size(v_map)
   nhist=sqrt(nx)
   iv_map=zeros(shape(v_map[lag:]))
   ivmax=abs(v_map.max())+abs(v_map.min())
   ivhist=array([0]*nhist)
   # find every relevant shifts and their number
   shift=lag
   inc=v_map-roll(v_map,shift)
   inc=inc[lag:]
   if is_angle:
      inc=within_pi2(inc)
   if (norm==2):
      iv_map=iv_map+inc**2
   else:
      iv_map=iv_map+abs(inc)
   ivmax=abs(inc).max()
   inc=array(list(inc)+list(-inc))
   if is_angle:
      nbs,bins=histogram(within_pi2(inc),bins=nhist,range=[-ivmax,ivmax])
   else:
      nbs,bins=histogram(inc,bins=nhist,range=[-ivmax,ivmax])
   ivhist=ivhist+array(nbs)
   cbin=array((bins[:-1]+bins[1:])*0.5)
   sig=sqrt(sum(ivhist*cbin**2)/sum(ivhist))
   if (norm==2):
      iv_map=sqrt(iv_map)
   return iv_map,cbin,double(ivhist),sig



# Saves array in nda format for DISPERSE (skeleton finding)
def save_nda(a,filename):
   from numpy import shape
   f=open(filename,'w')
   f.write('ANDFIELD\n')
   nx,ny=shape(a)
   f.write('['+str(nx)+' '+str(ny)+']\n')
   f.write('BBOX [-0.5,-0.5] [1.0,1.0]\n')
   for v in a.flatten():
      f.write(str(v)+'\n')
#   f.write(' '.join(["%.4f" % v for v in a.flatten()]))  #(str(a.flatten())[1:-1])
   f.close()


# Spline derivate a function
def derivate(v,t,k=4):
   from scipy import interpolate
   from numpy import zeros
   iv=interpolate.InterpolatedUnivariateSpline(t,v,k=k)
   h=1e-3*(t[2]-t[1])
   dv=zeros(len(v))
   dv[1:-1]=(iv(t[1:-1]+h)-iv(t[1:-1]-h))/(2.*h)
   dv[0]=(iv(t[0]+h)-v[0])/h
   dv[-1]=(v[-1]-iv(t[-1]-h))/h
   return dv


# Stop conditions
def stop_default(t,y):
   from pylab import norm
   return t>10. or norm(y)>10.

# Use dvode to solve ODEs dy/dt=f(t,y) with y=y0 at t=t0.
def dsolve(f,t0,y0,dt=0.05,stop_at=stop_default,ncalls_max=100000,verbose=False):
   from numpy import reshape,array
   from scipy.integrate import ode
   r=ode(f).set_integrator('vode',nsteps=100000)
   r.set_initial_value(y0,t0)
   ts=[t0]
   ys=[y0]

   ncalls=0
   while (r.successful() and (not stop_at(r.t,r.y)) and ncalls<ncalls_max):
      ncalls=ncalls+1
      r.integrate(r.t+dt)    # dt is the Sampling resolution in time
      ts.append(r.t)
      ys.append(r.y)
      if verbose:           
         print (ncalls,' y=',r.t)
   ys=reshape(ys,[len(ys),len(y0)])
   ts=array(ts)
   return ts,ys

