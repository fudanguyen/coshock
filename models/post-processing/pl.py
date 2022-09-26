#!/usr/bin/python
#
# This script loads and plots named 1D data
# Author: Pierre L.
#
# TODO:
# load a list of files.Done
# ..At a given position.
# merge two files into one list of data.
# Plot them with different styles
# markers
# filters
# recall plot from saved protocole ?
# Thick lines and big fonts
from scipy import *
from pylab import *
rcParams['axes.unicode_minus']=False
rc('lines', linewidth=3)
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 18}
rc('font', **font)
rcParams['legend.fontsize']='small'

# common data
ipos=0
npos=0
v=[{}]
dat_file='/tmp/pl.dat'
files=[]
line_styles=['-','--',':','-.','o-','o--','o:','o-.','^-','^--','^:','^-.','s-','s--','s:','s-.']
#line_styles=['-','--',':','-.']
colors='b g r c m y k'.split()
#colors=['k'] # for black and white

# save data
def save_commons():
    import pickle
    f=open(dat_file,'wb')
    pickle.dump(v,f)
    pickle.dump(ipos,f)
    pickle.dump(npos,f)
    pickle.dump(files,f)
    f.close()

# retrieve data
def retrieve_commons():
    global ipos,v,npos,files
    import pickle
    f=open(dat_file,'r')
    v=pickle.load(f)
    ipos=pickle.load(f)
    npos=pickle.load(f)
    files=pickle.load(f)
    f.close()
    
#----------------------
# file load
#----------------------

# Expand file_spec and read corresponding files
def fl(file_specs,mesa=False):
##     # retrieve data
##     import pickle
##     f=open(dat_file,'r')
##     v=pickle.load(f)
##     ipos=pickle.load(f)
##     npos=pickle.load(f)
##     files=pickle.load(f)
##     f.close()

    # parse file_spec
    import os
    for file_spec in file_specs.split(','):
        import re
        if re.search('//',file_spec):
            file_list=[file_spec]
        else:
            file_list=os.popen("ls "+file_spec).read().split()
        for file in file_list:
            file_load(file,mesa=mesa)

# saves a given position into a file
def file_save(filename,pos=-1):
    global ipos
    if (pos != -1 and pos < npos):
        ipos=pos
    f=open(filename,'w')
    print ('save file:',filename,' from pos ',ipos)
    keys=v[ipos].keys()
    print >>f,' '.join(keys)
    for i in range(len(v[ipos][keys[0]])):
        line=''
        for k in keys:
            line=line+str(v[ipos][k][i])+' '
        print>>f,line
    f.close()

# loads a file at a given position 
def file_load(filenames,mesa=False,pos=-1):
    global ipos,npos
    import mytools
    from numpy import array
    # Determines position where to load file
    if pos<0 :
        ipos=npos
        npos=npos+1
        v.append({})
        files.append(filenames)
    else:
        ipos=pos
        files[ipos]=filenames

        
    # Read data
    lines=[]
    for filename in filenames.split('//'):
        f=open(filename,'r')
        dlines=f.readlines()
        if lines==[]:
            lines=dlines
        else:
            # Glue the two files side by side
            n=max([len(dlines),len(lines)])
            m=min([len(dlines),len(lines)])-1
            if len(lines)>=len(dlines):                
                for i in range(n):
                    lines[i]=lines[i]+' '+dlines[min([i,m])]
            else:
                for i in range(n):
                    dlines[i]=lines[min([m,i])]+' '+dlines[i]
                lines=dlines
        f.close()
   
    # Skip nskip first lines (5 in MESA format)
    nskip=0
    if (mesa):
        nskip=5
    # loads data
    if options.no_header:
        nvars=len(lines[0].split())
        vars=map(lambda i:'v'+str(i),range(nvars))
        nskip=-1
    else:
        vars=lines[nskip].split()
    
    # Convert annoying characters
    for i in range(len(vars)):
        vv=vars[i]
        vv=vv.replace(',','_')
        vv=vv.replace(';','_')
#        vv=vv.replace('(','')
#        vv=vv.replace(')','')
        vars[i]=vv
#        print 'variable ',i,' is :',vv

    # verbose
    print ('load file:',filename,' at pos ',ipos)
#    print 'ipos:',ipos
#    print 'npos:',npos
 
    for var in vars:
        v[ipos][var]=[]
    import re
    for line in lines[1+nskip:]:
        vals=line.split()
        if len(vals)!=len(vars):
            print()
            print (line)
            print ('file ',filename)
            print ('This line has not the right amount of values:')
            print (len(vals),' values for ',len(vars),' variables.')
            for i in range(min([len(vals),len(vars)])):
                print ('var ',i,vars[i],' val=',vals[i])
            if (len(vals)>len(vars)):
                print ('extra values:')
                print (vals[len(vars):])
            if (len(vals)<len(vars)):
                print ('extra variables:')
                print (vars[len(vals):])
            import sys
            sys.exit()
        for (var,val) in zip(vars,vals):
            if re.search('[^E]-\d\d\d',val):
                val='0'
            if re.search('NaN',val):
             #   print 'I found a NaN for var:',var
                val='0'
            if mytools.type_of_value(val)==str:
            #    val='0'
                v[ipos][var].append(val)
            else:
                v[ipos][var].append(float(val))
			
    for k in v[ipos].keys():
        v[ipos][k]=array(v[ipos][k])
    save_commons()
    return ipos

#----------------------
# tools
#----------------------

# Expand the names in formula to get nice names
def nice_expand(form):
    import re
    if re.compile('[<>]').search(form):
        words=re.findall("<([^>]*)>",form)
        words=list(set(words))
        for w in words:
            fvar= expand(w.strip())
            form=re.sub('<'+w+'>',fvar,form)
        return form
    else:
        return expand(form)
    
# Extract vars from formula and replace them by  v[j][key][i] 
def form_expand(form):
    import re
    if re.compile('[<>]').search(form):
        # Add numpy. to functions
#        words=re.findall("[^<][^>]*",form)
        words=re.findall("[^<]([a-z]+)\(",form)
        words=list(set(words))
        for w in words:
            form=re.sub(w+'\(','numpy.'+w+'(',form)
        # expand <variable names>
        words=re.findall("<([^>]*)>",form)
        words=list(set(words))
        for w in words:
            fvar='v[j]["'+expand(w.strip())+'"][i]'
            form=re.sub('<'+w+'>',fvar,form)
        return form
    else:
        return 'v[j]["'+expand(form)+'"][i]'

# Evaluate formula on all indices of array, returns the array
def form_eval(form):
    import numpy 
    # Get size of array
    k=v[ipos].keys()
    n=len(v[ipos][k[0]])
    # Defines and fill the array
    arr=[0]*n
    j=ipos
    for i in range(n):
        arr[i]=eval(form)
    return numpy.array(arr)

# Find key corresponding to var
def expand(var):
    # If var is as such in keys list
    if var in v[ipos].keys():
        return var
    # Expand variable name
    import re
    short=re.compile(var)
    for k in v[ipos].keys():
        if short.match(k):
            return k
    print ("\n************\n variable ",var," not found in ")
    print (v[ipos].keys())
    print ('file:',files[ipos],"\n************\n")
        
#----------------------
# multi-plot
#----------------------
def mp(string):
    import pylab as p
    global ipos,npos,options
    
    retrieve_commons()
    
    # open figure
    p.figure(options.figure)

    # Thick lines
    p.rc('lines', linewidth=2)

    # Define log axis
    if (options.yl):
        p.semilogy()
    if (options.xl):
        p.semilogx()
    if (options.xl and options.yl):
        p.loglog()
        
    # For each file
    for image in range(npos):
        ipos=image
        if options.verbose:
            print ('file:',files[ipos])
        # parse ord;abs
        (ord,abs)=string.split(';')

        # define titles
        # defaults
        stitle=options.filename
        ytitle=ord
        xtitle=nice_expand(abs)
        # overwrite with options
        if (options.title):
            stitle=options.title
        if (options.xtitle):
            xtitle=options.xtitle
        if (options.ytitle):
            ytitle=options.ytitle
        p.xlabel(xtitle)
        p.ylabel(ytitle)
        #p.title(stitle)

        # parse ordinates 
        abs=form_expand(abs)
        ord=ord.split(',')
        nice_ord=map(nice_expand,ord)
        ord=map(form_expand,ord)
        if options.verbose:
            print (abs)

        # plot
        for (o,lo) in zip(ord,nice_ord):
            if options.verbose:
                print (o)
            if options.color:
                line_style=line_styles[ord.index(o) % len(line_styles)]
                color=colors[ipos % len(colors)]
            else:
                line_style=line_styles[ipos % len(line_styles)]
                color=colors[ord.index(o) % len(colors)]
            # Info on file except if too big
            lf=''
            label=None
            if (len(files[ipos])<15):
                lf=' ('+files[ipos]+')'
                label=lo+lf
            else:
                # If you don't info on file, then label only vars.
                if (ipos==0):
                    label=lo
            p.plot(form_eval(abs),form_eval(o),line_style,label=label,color=color)
    # Define and set ranges
    if (options.x):
        xrange=map(float,options.x.split(','))
        p.xlim(xrange)
    if (options.y):
        yrange=map(float,options.y.split(','))
        p.ylim(yrange)        

    # legend
    p.legend(loc=0)
    p.figure(options.figure).subplots_adjust(bottom=0.15)
    p.show()
    
# Defines and parse options, call file load or multi plot accordingly
from optparse import OptionParser
parser = OptionParser()
# Defines arguments
parser.add_option("-f", "--file",
                  dest="filename",
                  help="filename with data")
parser.add_option("--nh",
                  dest="no_header", default=False,action='store_true',
                  help="No header with the variable names: use v1, v2 ..")
parser.add_option("-p",
                  dest="plot_only", default=False,action='store_true',
                  help="plot only mode")    
parser.add_option("-l",
                  dest="load_only", default=False,action='store_true',
                  help="load only mode")    
parser.add_option("-m",
                  dest="mesa", default=False,action='store_true',
                  help="data file is MESA format")    
parser.add_option("-v",
                  dest="verbose", default=False,action='store_true',
                  help="toggles verbose mode")    
parser.add_option("-c",
                  dest="color", default=False,action='store_true',
                  help="toggles colors lead line_style")    
parser.add_option("-n",
                  dest="figure", default=1,
                  help="define figure number")    
parser.add_option("-X",
                  dest="xl", default=False,action='store_true',
                  help="log x axis")    
parser.add_option("-Y",
                  dest="yl", default=False,action='store_true',
                  help="log y axis")
parser.add_option("-x",
                  dest="x", default=False,
                  help="x axis range")    
parser.add_option("-y",
                  dest="y", default=False,
                  help="y axis range")
parser.add_option("-t","--title",
                  dest="title", default=False,
                  help="figure title")    
parser.add_option("--xt","--xtitle",
                  dest="xtitle", default=False,
                  help="figure x title")    
parser.add_option("--yt","--ytitle",
                  dest="ytitle", default=False,
                  help="figure y title")    
# Default options
(options,args)=parser.parse_args([])

# Take options from command line arguments
if __name__ == "__main__":
    # Parse them
    (options, args) = parser.parse_args()

    # verbose
    if (options.verbose):
        print ('filename:',options.filename)
        print ('graph request:',args[0])
        print ('figure number:',options.figure)

    # actions
    if (options.load_only or not options.plot_only):
        fl(options.filename,mesa=options.mesa)
    if (options.plot_only or not options.load_only):
        mp(args[0])

	
	

