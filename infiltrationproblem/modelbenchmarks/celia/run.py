import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from celia_MPM import ModelRun
import time

# MyNumba.py
from numba import types
from numba.typed import Dict

def MakeDictFloat():
    d=Dict.empty(
    key_type=types.unicode_type,
    value_type=types.float64,)
    return d

# Driving data:
qT=np.loadtxt('../../input/infiltration.dat',skiprows=1,delimiter=',',usecols=1)/1000*-1.
qB=np.zeros(len(qT))
t=np.arange(len(qT))
dt=t[1]-t[0]
nt=len(t)

#qT=qT[:1000]
#t=t[:1000]
#nt=len(t)

# Soil properties:
pars=MakeDictFloat()
#pars={}
pars['thetaR']=0.131
pars['thetaS']=0.396
pars['alpha']=0.423
pars['n']=2.06
pars['m']=1-1/pars['n']
pars['Ks']=0.0496
pars['neta']=0.5
pars['Ss']=0.000001

# Spatial grid:
dz=0.1
zN=1.5
z=np.arange(dz/2,zN,dz)
n=len(z)

# Initial condition:
psi=np.zeros((nt,n))-3.59

psiB=np.array([-3.59])
psiT=np.array([-3.59])


# Run model m times
for repeat in range(2):
    tic=time.time()
    psi,QIN,QOUT,S,err=ModelRun(dt,dz,n,nt,psi,qT,psiB,psiT,pars)
    runtime=time.time()-tic
    print('Runtime = %.4f seconds'%(runtime))

# Pack output into a dataframe:
WB=pd.DataFrame(index=t)
WB['S']=S
WB['QIN']=-QOUT
WB['QOUT']=-QIN


psiRMSE=0.0

print(psi.shape)
print(WB.shape)
np.save('psi.npy',psi)

# Process and save output
t=WB.index
S=WB['S'].values
Qin=WB['QIN'].values
Qout=WB['QOUT'].values

dt=t[1]-t[0]
dS=S[-1]-S[0]
dQ=np.sum(Qin-Qout)*dt

QIN=np.sum(Qin)*dt*1000
QOUT=np.sum(Qout)*dt*1000
err=1-dQ/dS
MBerr=(dS-dQ)*1000
RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt)**2))*1000
MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt))*1000

np.save('WB',WB.values)

f=open('result.csv','w')
f.write('Q approx, Runtime (s), psiRMSE (m), MB err (mm), MB RMSE (mm), QIN (mm), QOUT (mm), dS (mm)\n')
f.write('central, %.3f, %.3f, %.3f, %.3e, %.2f %.2f, %.2f\n'%(runtime, psiRMSE, MBerr, RMSE, QIN, QOUT, dS*1000))
try: 
    Qin=WB['QIN_FD'].values
    Qout=WB['QOUT_FD'].values
    QIN=np.sum(Qin)*dt*1000
    QOUT=np.sum(Qout)*dt*1000
    err=1-dQ/dS
    MBerr=(dS-dQ)*1000
    RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt)**2))*1000
    MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt))*1000
    f.write('forward, %.3f, %.3f, %.3f, %.3e, %.2f %.2f, %.2f\n'%(runtime, psiRMSE, MBerr, RMSE, QIN, QOUT, dS*1000))
except:
    # Nothing to do
    print('OK') 

try:
    Qin=WB['QIN_BD'].values
    Qout=WB['QOUT_BD'].values
    QIN=np.sum(Qin)*dt*1000
    QOUT=np.sum(Qout)*dt*1000
    err=1-dQ/dS
    MBerr=(dS-dQ)*1000
    RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt)**2))*1000
    MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt))*1000
    f.write('backward, %.3f, %.3f, %.3f, %.3e, %.2f %.2f, %.2f\n'%(runtime, psiRMSE, MBerr, RMSE, QIN, QOUT, dS*1000))
except:
    # Nothing to do
    print('OK') 

f.close()
