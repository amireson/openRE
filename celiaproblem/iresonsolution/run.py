import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from Richards import run_RE
from Richards import KFun

# Run the models:
def KFun(psi,pars):
    x2=pars['A']+np.abs(psi)**pars['gamma']
    x1=pars['A']/x2
    K=pars['Ks']*x1
    return K


# Soil properties:
def setpars():
    pars={}
    pars['thetaR']=0.075
    pars['thetaS']=0.287
    pars['alpha']=1.611e6
    pars['beta']=3.96
    pars['A']=1.175e6
    pars['gamma']=4.74
    pars['Ks']=0.00944
    return pars

pars=setpars()

# Time grid:
tN=360.

# Spatial grid:
dz=1.
zN=40.
z=np.arange(dz,zN,dz)
n=len(z)
z=np.hstack([0,z,zN])
#z=z[-1]-z

# Initial condition:
psi0=np.zeros(n)-61.5

# Boundary conditions:
psiT=-20.7
psiB=-61.5

dti=np.array([0.1,1,3,10,20,30,40,60,90,120])
for dt in dti:
    t=np.arange(0,tN+dt,dt)
    nt=len(t)
    BC_T=psiT+np.zeros(nt)
    BC_B=psiB+np.zeros(nt)
    psi,WB,runtime=run_RE(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)

    zs=np.zeros(nt)
    Kin=(KFun(np.array([psiT]),pars)+KFun(psi[:,0],pars))/2.
    qin=-Kin*((psi[:,0]-psiT)/dz-1.)
    Kout=(KFun(psi[:,-1],pars)+KFun(np.array([psiB]),pars))/2.
    qout=-Kout*((psiB-psi[:,-1])/dz-1.)

    WB['QIN_BD']=zs
    WB['QIN_BD'].iloc[1:]=qin[1:]

    # Lower BC: Type 1 specified pressure head
    WB['QOUT_BD']=zs
    WB['QOUT_BD'].iloc[1:]=qout[1:]

    psi=np.hstack([psiT,psi[-1,:],psiB])

    # Save data:
    dtstr=str(dt).replace('.','-')
    fname='output/psi_dt%s.csv'%dtstr
    f=open(fname,'w')
    f.write('z, psi\n')

    for i,j in zip (z,psi): f.write('%.4f, %.4f\n'%(i,j))
    f.close()
 
    dtstr=str(dt).replace('.','-')
    fname='output/mb_dt%s.csv'%dtstr
    WB.to_csv(fname,float_format='%12.8f')
