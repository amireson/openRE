import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from celia_PI import ModelRun

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

def setup(dt):
    # Set parameters:
    pars=setpars()

    # Grid:
    zN=40.
    dz=1.
    tN=360.

    z=np.arange(dz,zN,dz)
    n=len(z)

    t=np.arange(0,tN+dti,dti)
    nt=len(t)

    # Initialize array:
    psi=np.zeros((nt,n))

    # ICs:
    psi[0,:]=-61.5

    # BCs:
    psiB=-61.5
    psiT=-20.7

    return z,t,dz,n,nt,zN,psi,psiB,psiT,pars

dt=np.array([0.1,1,3,10,20,30,40,60,90,120])
for dti in dt:
    z,t,dz,n,nt,zN,psi,psiB,psiT,pars=setup(dti)
    psi,QIN,QOUT,S,err=ModelRun(dti,dz,n,nt,psi,psiB,psiT,pars)
    z=np.hstack([0,z,zN])
    z=z[-1]-z
    psi=np.hstack([psiB,psi[-1,:],psiT])

    # Save data:
    dtstr=str(dti).replace('.','-')
    fname='output/psi_dt%s.csv'%dtstr
    f=open(fname,'w')
    f.write('z, psi\n')
    for i,j in zip (z,psi): f.write('%.4f, %.4f\n'%(i,j))
    f.close()
    
    dtstr=str(dti).replace('.','-')
    fname='output/mb_dt%s.csv'%dtstr
    WB=pd.DataFrame(index=t)
    WB['S']=S
    WB['QIN_FD']=np.nan
    WB['QOUT_FD']=np.nan
    # These are reversed, since I am using z as depth, Celia uses elevation
    WB['QIN_BD']=-QOUT
    WB['QOUT_BD']=-QIN
    WB['QIN_CN']=np.nan
    WB['QOUT_CN']=np.nan
    WB.to_csv(fname,float_format='%12.8f')
