import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import time

import sys
sys.path.insert(0,'..')
from celia_MPM_T1_VG import ModelRun

from numba import types
from numba.typed import Dict

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

from setup import *
pars_dict=pars
pars=MakeDictFloat()
for key in pars_dict: pars[key]=pars_dict[key]

z=np.arange(dz,zN,dz)
n=len(z)

# Initial condition:
psi0=np.zeros(n)+psiB

t=np.arange(0,tN+dt,dt)
nt=len(t)
#BC_T=psiT+np.zeros(nt)
psiT=np.array([psiT])
psiB=np.array([psiB])

to=[5,10,20,50,100]
tint=np.rint(t*24*60).astype('int')
ti=[np.where(tint==i)[0][0] for i in to]

#BC_B=psiB+np.zeros(nt)

psi=np.zeros((nt,n))
psi[0,:]=psi0

tic=time.time()
dum,QIN,QOUT,S,err=ModelRun(dt,dz,n,2,psi,psiB,psiB,pars)
runtime=time.time()-tic
print(runtime)

tic=time.time()
psi,QIN,QOUT,S,err=ModelRun(dt,dz,n,nt,psi,psiT,psiB,pars)
runtime=time.time()-tic
print(runtime)


theta=thetaFun(psi[ti,:],pars)
np.savetxt('theta.csv', np.vstack([z,theta]).T, delimiter=',',fmt='%10.4f')

psi=np.hstack([psiT,psi[-1,:],psiB])

# Save data:
f=open('runtime.csv','w')
f.write('%.3f\n'%runtime)
f.close()

WB=pd.DataFrame(index=t)
WB['S']=S
WB['QIN']=QIN
WB['QOUT']=QOUT

fname='mb.csv'
WB.to_csv(fname,float_format='%12.8f')
