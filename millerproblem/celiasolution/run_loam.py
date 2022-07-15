import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from numba import jit 
import time

#import sys
#sys.path.insert(1,'../../library/celia')
from celia_MPM_T1_VG_psiall import ModelRun
#from celia_MPM_original import ModelRun

# MyNumba.py
from numba import types
from numba.typed import Dict

dt=0.001

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

def setpars():
    pars=MakeDictFloat()
    pars['thetaR']=0.078
    pars['thetaS']=0.430
    pars['alpha']=3.600
    pars['n']=1.560
    pars['m']=1-1/pars['n']
    pars['Ks']=0.250
    pars['neta']=0.5
    pars['Ss']=1e-6
    return pars

def setup(dt):
    # Set parameters:
    pars=setpars()

    # Grid:
    dz=0.0125
    zN=5.

    tN=2.25

    z=np.arange(dz,zN,dz)
    n=len(z)

    t=np.arange(0,tN+dt,dt)
    nt=len(t)

    # Initialize array:
    psi=np.zeros((nt,n))

    # ICs:
    psi[0,:]=-z

    # BCs:
    psiT=np.array([0.1])
    psiB=np.array([0])

    return z,t,dz,n,nt,zN,psi,psiB,psiT,pars,dt

z,t,dz,n,nt,zN,psi,psiB,psiT,pars,dt=setup(dt)

tic=time.time()
dum,QIN,QOUT,S,err=ModelRun(dt,dz,n,2,psi,psiB,psiT,pars)
runtime=time.time()-tic
print('Celia solution, dt=%.4f, runtime = %.2f seconds'%(dt,runtime))

tic=time.time()
psi,QIN,QOUT,S,err=ModelRun(dt,dz,n,nt,psi,psiB,psiT,pars)
runtime=time.time()-tic
print('Celia solution, dt=%.4f, runtime = %.2f seconds'%(dt,runtime))

z=np.hstack([0,z,zN])
z=z[-1]-z
psiB=np.zeros((nt,1))+psiB
psiT=np.zeros((nt,1))+psiT
psi=np.hstack([psiB,psi,psiT])

# Save data:
dtstr=str(dt).replace('.','-')
fname='output_loam/psi_dt%s.csv'%dtstr
f=open(fname,'w')
f.write('z, psi\n')
for i,j in zip (z,psi[-1,:]): f.write('%.4f, %.4f\n'%(i,j))
f.close()

dtstr=str(dt).replace('.','-')
f=open('output_loam/runtime_dt%s.csv'%dtstr,'w')
f.write('%.3f'%runtime)
f.close()

fname='output_loam/mb_dt%s.csv'%dtstr
f=open(fname,'w')
f.write('t,S,QIN,QOUT\n')
for i,j,k,l in zip (t,S,QIN,QOUT): f.write('%.4f, %.8f, %.8f, %.8f\n'%(i,j,k,l))
f.close()

