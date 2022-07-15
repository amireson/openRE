import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

dt='0-001'

def GetError(df):
    dt=df.index[1]-df.index[0]
    Qin=df['QIN'].values
    Qout=df['QOUT'].values
    S=df['S'].values
    
    dS=S[-1]-S[0]
    dQ=np.sum((Qin-Qout)*dt)
    
    QIN=np.sum(Qin)*dt
    QOUT=np.sum(Qout)*dt
    err=1-dQ/dS
    
    MBerr=(dS-dQ)
    RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt)**2))
    MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt))
    return MBerr, RMSE, QIN

pl.figure(figsize=(4.5,4))

#fname='celiasolution/CeliaMillerProb/output_sand/psi_dt%s.csv'%dti
#z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
#i=np.arange(0,len(z),10)
#pl.plot(psi[i],10-z[i],'--',color='grey')

fname='../iresonsolution/sand/output/psi.csv'
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
pl.plot(psi,10-z,'-',color='tab:blue',label='sand (rtol=1e-6, atol=1e-6)')

fname='../celiasolution/output_sand/psi_dt%s.csv'%dt
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
pl.plot(psi,10-z,'--',color='silver')

#fname='celiasolution/CeliaMillerProb/output_loam/psi_dt%s.csv'%dti
#z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
#i=np.arange(0,len(z),10)
#pl.plot(psi[i],5-z[i],'--',color='grey')

fname='../iresonsolution/loam/output/psi.csv'
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
pl.plot(psi,5-z,color='tab:orange',label='loam (rtol=1e-6, atol=1e-6)')

fname='../celiasolution/output_loam/psi_dt%s.csv'%dt
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
pl.plot(psi,5-z,'--',color='silver')

#fname='celiasolution/CeliaMillerProb/output_clayloam/psi_dt%s.csv'%dti
#z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
#i=np.arange(0,len(z),10)
#pl.plot(psi[i],2-z[i],'--',color='grey')

fname='../iresonsolution/clayloam/output/psi.csv'
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
pl.plot(psi,2-z,color='tab:green',label='clayloam (rtol=1e-5, atol=1e-5)')

fname='../celiasolution/output_clayloam/psi_dt%s.csv'%dt
z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
i=np.arange(0,len(z),10)
pl.plot(psi,2-z,'--',color='silver')

#pl.plot(np.nan,np.nan,'--',color='grey',label='Celia solutions, dt=0.01d')
pl.plot(np.nan,np.nan,'--',color='silver',label='Celia solutions, dt=0.001d')

pl.grid()
pl.xlabel('Pressure head (cm)')
pl.ylabel('Elevation (cm)')
pl.legend()
pl.subplots_adjust(top=0.96)
#pl.show()
pl.savefig('Figure04.png',dpi=300)
