#!/usr/bin/python
#
# strobe-X_simulate.py
#
# Simulate data using the fakeit command
# Uses the Strobe-X response files, and the relxill model
# 
# This version uses parameters appropriate for AGN 
#
# Requires: xspec
#
import sys 
from xspec import *
from optparse import OptionParser
import os,os.path
import glob
import numpy as np
from subprocess import call
#
# ------------------------------------------------------------------------------
#
# MAIN PROGRAM
#
#
#
version='0.1a'
date='- Tue Sep 12 16:06:53 PDT 2017 -'
author='Javier Garcia'
#
ul=[]
ul.append("usage: %prog [options] PREFIX")
ul.append("")
ul.append("Get total counts in different bands for a given observation")
ul.append("PREFIX can be a single PHA file or a group (e.g. *.pha)")
usage=""
for u in ul: usage+=u+'\n'

parser=OptionParser(usage=usage)
parser.add_option("-v","--version",action="store_true",dest="version",default=False,help="show version number")
parser.add_option("-i","--input",dest="inpfile",default="",help="specify spectrum file (default: meg_-1.pha.gz)")
parser.add_option("-o","--output",dest="outfile",default="",help="specify alternative output file (default: counts.out)")

(options,args)=parser.parse_args()

#-----
# No chatter
Xset.chatter = 0

# Query
Fit.query = 'yes'

# Set abundances and cross sections
Xset.abund = "wilm"
Xset.xsect = "vern"

# Response files
respath="/Users/javier/simulations/StrobeX/responses/"
res1="xrca_2017-08-11.rmf"
arf1="xrca_128_2017-08-11.arf"
back1=""
res2="LAD_90mod_200eV.rmf"
arf2="LAD_90mod_200eV.arf"
back2="LAD_90mod_200eV.bkg"

# Create symbolic links
call(["ln","-s",respath+res1])
call(["ln","-s",respath+arf1])
call(["ln","-s",respath+res2])
call(["ln","-s",respath+arf2])
call(["ln","-s",respath+back2])

# Parameter of interest
kTe_vals=['20','40','60','80','100','120','140','160','180','200']

# Numer of simulations per case
maxiter = 10

# Source and Background exposure times
stime=2.e4
btime=stime

# Source flux
f1mC = 3.e-11  # Flux for 1 mCrab
sflux = 10.*f1mC

# Load local models
AllModels.lmod("relxill")

# Output
outfile='spfit-kTe.error'
f = open(outfile, 'w')

for kTe in kTe_vals:


  for it in range(maxiter):

    print "*** kTe",kTe,"iter",it

    # Define the Model
    m1 = Model("tbabs*relxillCp")

    m1(1).values = "0.5 -1"      # Tbabs   Nh
    m1(2).values = "5. -1"       # relxill Index1
    m1(3).values = "3. -1"       # relxill Index2
    m1(4).values = "5. -1"       # relxill Rbr
    m1(5).values = "0.998 -1"    # relxill a
    m1(6).values = "45. 1"       # relxill Incl
    m1(7).values = "-1. -1"      # relxill Rin
    m1(8).values = "400. -1"     # relxill Rout
    m1(9).values = "1.e-2 -1"    # relxill z
    m1(10).values = "2.0 1"      # relxill gamma
    m1(11).values = "2. 1"       # relxill logxi
    m1(12).values = "3. -1"      # relxill Afe
    m1(13).values = kTe+" 1"     # relxill kTe
    m1(14).values = "1.0 -1"     # relxill refl_frac
    m1(15).values = "1. 1"       # relxill norm

    # Calculate 2-10 keV flux
    AllModels.calcFlux("2. 10.")
    flux=AllModels(1).flux[0]

    # New normalization for 1 mCrab
    norm = sflux/flux
    m1(15).values = str(norm)      # relxill norm

    # Reset variables
    cval=[]
    wval=[]

    # First simulate XRCA 
    #specfile
    specfile='sim_xrca_128_kTe-'+kTe+'.fak'
    grpspecfile='sim_xrca_128_kTe-'+kTe+'_grp.fak'
    #response, arf, background, exposure, correction, backExposure, fileName
    fs1 = FakeitSettings(res1,arf1,back1,stime,1.,btime,specfile)
    AllData.fakeit(1, fs1)

    # Unload data
    AllData.clear()

    # Then simulate LAD
    #specfile
    specfile2='sim_LAD_90mod_200eV_kTe-'+kTe+'.fak'
    grpspecfile2='sim_LAD_90mod_200eV_kTe-'+kTe+'_grp.fak'
    #response, arf, background, exposure, correction, backExposure, fileName
    fs2 = FakeitSettings(res2,arf2,back2,stime,1.,btime,specfile2)
    AllData.fakeit(1, fs2)

    # Unload data
    AllData.clear()

    # Rebin the data
    call(["cp",specfile,grpspecfile])
    call(["/Users/javier/bin/isis","rebin_spec_simple.sl",grpspecfile])

    call(["cp",specfile2,grpspecfile2])
    call(["/Users/javier/bin/isis","rebin_spec_simple.sl",grpspecfile2])

    # Unload the model
    AllModels.clear()

    #
    # Load data
    AllData('1:1 '+grpspecfile+' 2:2 '+grpspecfile2)

    s1 = AllData(1)
    s2 = AllData(2)

    # Ignore/notice data
    s1.ignore("0.-0.2,12.-**")
    s2.ignore("0.-2.,30.-**")

    # Define the Model
    AllModels += "tbabs*relxillCp"
    m1 = AllModels(1)
    m2 = AllModels(2)

    m1(1).values = "0.5 -1"      # Tbabs   Nh
    m1(2).values = "5. -1"       # relxill Index1
    m1(3).values = "3. -1"       # relxill Index2
    m1(4).values = "5. -1"       # relxill Rbr
    m1(5).values = "0.998 -1"    # relxill a
    m1(6).values = "45. 1"       # relxill Incl
    m1(7).values = "-1. -1"      # relxill Rin
    m1(8).values = "400. -1"     # relxill Rout
    m1(9).values = "1.e-2 -1"    # relxill z
    m1(10).values = "2.0 1"      # relxill gamma
    m1(11).values = "2. 1"       # relxill logxi
    m1(12).values = "3. 1"       # relxill Afe
    m1(13).values = kTe+" 1"     # relxill kTe
    m1(14).values = "1.0 1"      # relxill refl_frac
    m1(15).values = "1. 1"       # relxill norm

    m2(15).values = "=15"

    # Fit
    Fit.renorm()
    Fit.perform()

    # Calculate Errors
    Fit.error("maximum 3. 2.706 13")

    # Save values
    cval.append(m1(13).values[0])
    wval.append(1./((m1(13).error[1]-m1(13).error[0])/2.)**2.) # There should be a better way!

    # Unload data
    AllData.clear()

    # Unload the model
    AllModels.clear()

  # Calculate the weighted average and error
  avgVal = np.average(np.array(cval),weights=np.array(wval),returned=True)

  # Output
  f.write(str(kTe)+' '+str(avgVal[0])+' '+str(np.sqrt(1./avgVal[1]))+'\n')

f.close()
#
#
sys.exit()
# ------------------------------------------------------------------------------
