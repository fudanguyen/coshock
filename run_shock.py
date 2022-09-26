# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Mon Apr 12 00:31:34 2021

# @author: datng
# """

# =============================================================================
# Shell script to run shock models
# =============================================================================

import os
# https://stackoverflow.com/questions/431684/equivalent-of-shell-cd-command-to-change-the-working-directory/13197763#13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
        
import sys
import subprocess
import pexpect

shockRoutine = True
coRoutine = True
# Should be both True for new models run
# False-True for rerun of CO routine

# INPUT SET1
path_orig = '/home/datng/Documents/ic443/reach_2019/model_3_2021'
path_shock = path_orig+'/Shock_1.2'
path_input = '/input/input_mhd.in'
path_species = '/input/species.in'

shock_type = ['C']
# PHYSICAL INPUT
Bbeta = ['3.50E+00']
Vs_km = ['6.00E+01'] #['4.00E+01', '3.50E+01','4.50E+01']
nH_init = ['2.00E+03']
Tn = '2.00E+01'
RAD = '1.0'
# RUNTIME INPUT
Nstep_max = '18000'
timeJ = ['5E+03']
duration_max = '1E+06'

strings_param = [
    '! modele', '! specfile', '! chemfile',
    '! shock_type', '! Nfluids',
    '! Bbeta', '! Vs_km', '! nH_init', '! Tn',
    '! RAD', '! Cool_KN',
    '! Nstep_max', '! timeJ', '! duration_max']

with open(path_orig+'/input_mhd_blank.in', 'r') as input_file:
    strings_id, skips = [], []
    param = {}
    for i, line in enumerate(input_file):
        for match in strings_param:
            if (match in line) and not (match in skips):
                strings_id.append(i)
                skips.append(match)
                param[match] = i
                break
designated_spaces = 41

# # create input file for each shock model:
# open reference input file
# for type in shock_type:
# 	for n in nH_init:
# 		for b in Bbeta:
# 			for Vs_km in Vs_km:
# 				for Tn in Tn:
# 					for rad in RAD:
# for y in timeJ:
sfn_arr=[]
if True:
    for typ in shock_type:
    # 
        # for n in nH_init:
        if True:
            # for V in Vs_km:
            for y,n,V in zip(timeJ,nH_init,Vs_km): 
            
                #### Static model
                str_nH = "{0:1.0E}".format(int(float(n)))
                # str_b = "{0:1.1F}".format(float(b))
                Sname1 = 'S1'+'n'+str_nH[0]+'e'+str_nH[-1]
                # Sname2 = 'y'+str(y[0])+str(y[-1])+'v'+str(int(float(V)))
                Sname2 = 'y'+str(y[-1])+'v'+str(int(float(V)))
                Sfilename = path_orig+"/input_files/S/"+Sname1+Sname2
                sfn_arr.append(Sfilename)
                # Copy reference file to new input file
                if shockRoutine:
                    subprocess.call(["cp", path_orig+"/input_mhd_blank.in", 
                         Sfilename])
                
                listS = [
                    Sname1+'/'+Sname2, 'species.in_depl', 'chemistry.in_noadso',
                    'S1', '1',
                    '0.00E+00', V, n, Tn, 
                    RAD, '0', 
                    Nstep_max, y, duration_max]
                
                # Create dict. for writing values.
                Sinput = {}
                for i, names in enumerate(strings_param):
                    Sinput[names] = listS[i]
                
                # Write target file
                hdl = strings_param.copy()
                if shockRoutine:
                    with open(Sfilename, 'r') as Sfile:
                        lines = Sfile.readlines()
                        
                        for names in strings_param:
                            inputval = Sinput[names]
                            lines[param[names]] = inputval+(41-len(inputval))*' '+lines[param[names]]
                if shockRoutine:
                    with open (Sfilename, 'w+') as Sfile:
                        Sfile.writelines([line for line in lines])
                    
                    # Copy to new input file, run static model
                    subprocess.call(["cp", Sfilename, path_shock+path_input])
                    with cd('Shock_1.2'):
                        subprocess.call(["mkdir","output/"+Sname1+'/'+Sname2,"-p"])
                        subprocess.call(["./mhd_vode"])
                        
                    subprocess.call(["cp",path_shock+"/output/"+Sname1+'/'+Sname2+"/species.out",
                                         path_shock+"/input/species.in"])
                
                for b in Bbeta:
                    str_b = "{0:1.1F}".format(float(b))
                    Mname1 = typ+'n'+str_nH[0]+'e'+str_nH[-1]+'b'+str_b
                    Mname2 = 'y'+str(y[-1])+'v'+str(int(float(V)))
                    Mfilename = path_orig+"/input_files/"+Mname1+Mname2
                    
                    if shockRoutine:
                        # Copy reference file to new input file
                        subprocess.call(["cp", path_orig+"/input_mhd_blank.in", 
                                 Mfilename])
                        
                    if typ == "J": Nflu = '1'
                    else: Nflu = '3'
                    listM = [
                        Mname1+'/'+Mname2, 'species.in', 'chemistry.in_standard',
                        typ, Nflu,
                        b, V, n, Tn, 
                        RAD, '1', 
                        Nstep_max, y, duration_max]
                    
                    # Create dict. for writing values.
                    Minput = {}
                    for i, names in enumerate(strings_param):
                        Minput[names] = listM[i]
                    
                    # Write target file
                    if shockRoutine:
                        with open(Mfilename, 'r') as Mfile:
                            lines = Mfile.readlines()
                            for names in strings_param:
                                inputval = Minput[names]
                                lines[param[names]] = inputval+(41-len(inputval))*' '+lines[param[names]]
                        with open (Mfilename, 'w+') as Mfile:
                            Mfile.writelines([line for line in lines])
                        
                        # Copy to new input file, run static model
                        subprocess.call(["cp", Mfilename, path_shock+path_input])
                        with cd('Shock_1.2'):
                            subprocess.call(["mkdir","output/"+Mname1+'/'+Mname2,"-p"])
                            subprocess.call(["./mhd_vode"])
                        
                    if coRoutine:
                        child = pexpect.spawn('./../line_profiler/CO_cgrid-120619.exe')
                        child.expect('Then write : n4-b1')
                        child.sendline(Mname1+'\n')
                        child.expect('write : v30')
                        child.sendline(Mname2+'\n')
                        child.logfile=sys.stdout
                        child.expect(pexpect.EOF, timeout=None)
                