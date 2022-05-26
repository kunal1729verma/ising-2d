from __future__ import print_function

from mpi4py import MPI
import copy
import numpy as np
from scipy.special import comb
from numpy import linalg as la
import json
from pprint import pprint
import sys
import os
import os.path


def int_to_string(num,digs):
    return((str(num)).zfill(digs))

def float_to_string(num,digs,dec):
    tmpstring = ""
    fmtstring = "%"+"."+str(dec)+"f"
    tmpstring=fmtstring%num
    return((str(tmpstring).zfill(digs+dec+1)))


#--------------------------------------------------------------------
junktab="junk.run.junk"
runcasename="seql"

casename = "test"
L = 64
Tmin = 2.0e0
Tmax = 2.5e0
dT = 0.01e0
neqsweeps = 50000
nsamsweeps = 500000
nsampstp = 10
nblen = 200
initialize_all_up = False
print_cor = True
sequential = True
tmax = 500
Llist = [L]

#initialize MPI
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
nprocs = comm.Get_size()
MASTER = 0

print(len(sys.argv))
print(sys.argv)

if(myrank == 0) :
    #------------------------read input if it is there-----------------
    if(len(sys.argv) > 1) :
        runpfile = sys.argv[1]
        with open(runpfile) as runp_file:
            runpdata = json.load(runp_file)
            pprint(runpdata)
            for key in runpdata:
                if (key == 'casename') :
                    casename = runpdata['casename']
                elif (key == 'Llist') :
                    Llist = runpdata['Llist']
                elif (key == 'Tmin') :
                    Tmin = runpdata['Tmin']
                elif (key == 'Tmax') :
                    Tmax = runpdata['Tmax']
                elif (key == 'dT') :
                    dT = runpdata['dT']
                elif (key == 'neqsweeps') :
                    neqsweeps = runpdata['neqsweeps']
                elif (key == 'nsamsweeps') :
                    nsamsweeps = runpdata['nsamsweeps']
                elif (key == 'nsampstp') :
                    nsampstp = runpdata['nsampstp']
                elif (key == 'nblen') :
                    nblen = runpdata['nblen']
                elif (key == 'initialize_all_up') :
                    initialize_all_up = runpdata['initialize_all_up']
                elif (key == 'print_cor') :
                    print_cor = runpdata['print_cor']
                elif (key == 'sequential') :
                    sequential = runpdata['sequential']
                elif (key == 'tmax') :
                    tmax = runpdata['tmax']
                else :
                    print("Unknown card ", key)
                    exit()
                    
    rundat = {
        'casename' : casename,
        'Llist' : Llist,
        'Tmin' : Tmin,
        'Tmax' : Tmax,
        'dT' : dT,
        'neqsweeps' : neqsweeps,
        'nsamsweeps' : nsamsweeps,
        'nsampstp' : nsampstp,
        'nblen' : nblen,
        'initialize_all_up' : initialize_all_up,
        'print_cor' : print_cor,
        'sequential' : sequential,
        'tmax' : tmax
    }
    isingdat = {
        'casename' : casename,
        'L' : L,
        'Tmin' : Tmin,
        'Tmax' : Tmax,
        'dT' : dT,
        'neqsweeps' : neqsweeps,
        'nsamsweeps' : nsamsweeps,
        'nsampstp' : nsampstp,
        'nblen' : nblen,
        'initialize_all_up' : initialize_all_up,
        'print_cor' : print_cor,
        'sequential' : sequential,
        'tmax' : tmax               
    }
    
else :
    rundat = None
    isingdat = None

runcard = comm.bcast(rundat,root=MASTER)
isingcard = comm.bcast(isingdat,root=MASTER)

print(myrank,isingcard)
print(myrank,runcard)

ndigs=4
ndecs=4
dot="."
junk="jnk"
infile= dot+junk+dot+"in"+int_to_string(myrank,4)+".json."+junk
outfile = dot+junk+dot+"out"+int_to_string(myrank,4)+".out."+junk

print(myrank,infile)
print(myrank,outfile)


for L in runcard['Llist'] :
    isingcard['L'] = L
    Tlist = []
    Tvallist =[]
    lcount = 0
    Tval = runcard['Tmin']
    while Tval <= runcard['Tmax'] + (runcard['dT']/2.e0):
        Tvallist.append(Tval)
        Tval   += runcard['dT']
        lcount +=1

    nTs = len(Tvallist)
    nload = nTs//nprocs

    Tdone = 0;
    for iproc in range(0,nprocs):
        for iload in range(0,nload):
            if(myrank == iproc) :
                Tlist.append(Tvallist[iproc*nload+iload])
            Tdone+=1
    for il in range(Tdone,nTs):
        if(myrank == il-Tdone):
            Tlist.append(Tvallist[il])
            
    print(myrank, len(Tlist), Tlist)
    comm.Barrier()
    #print(myrank, len(Tlist), Tlist)
    lcasename =  runcard['casename']+"_L"+int_to_string(L,ndigs)
    #remove any previous versions
    removestring = "rm -f "+dot+lcasename+"*out.plt"
    os.system(removestring)
    collatestring = "cat "+dot+lcasename+"*out.plt > "+lcasename+"_out.plt"

    print(myrank,lcasename)
    for Tval in Tlist:
        runcase=dot+lcasename+"_T"+float_to_string(Tval,ndigs,ndecs)
        isingcard['casename'] = runcase
        isingcard['Tmin'] = Tval
        isingcard['Tmax'] = Tval+runcard['dT']/1000.e0
        isingcard['dT'] = 10.e0*runcard['dT']
        with open(infile, 'w') as inputfile:
            json.dump(isingcard, inputfile)
        runcommand = "./ising < "+infile + " > " +outfile
        #teststr = "echo '"+runcommand+"'"
        os.system(runcommand)
        if myrank == MASTER :
            os.system(collatestring)
    comm.Barrier()
    if myrank == MASTER :
        os.system(collatestring)
    comm.Barrier()
    
                




