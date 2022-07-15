import numpy as np
#import xarray as xr
import pandas as pd
import matplotlib.pyplot as pl
import os

def GetError(df):
    dt=df['dt']
    Qin=df['Qin'].values
    Qout=df['Qout'].values
    S=df['S'].values
    
    dS=S[-1]-S[0]
    dQ=np.sum((Qin-Qout)*dt)
    
    QIN=np.sum(Qin)*dt
    QOUT=np.sum(Qout)*dt
    err=1-dQ/dS
    
    MBerr=(dS-dQ)
    RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt[1:])**2))
    MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt[1:]))
    return MBerr, RMSE

def GetIresonData(mypath):
    os.system('cat ../%s/result.csv'%mypath)
    fname='../%s/WB.npy'%mypath
    data=np.load(fname)
    t=pd.date_range(start='1999-10-1',freq='D',periods=len(data))
    WB=pd.DataFrame(index=t)
    S=data[:,0]*1000
    Qin=data[:,1]*1000
    Qout=data[:,2]*1000
    
    dt=1.
    df=pd.DataFrame(index=t)
    df['dt']=dt
    df['Qin']=Qin
    df['Qout']=Qout
    df['S']=S
    df['dS']=np.hstack([0,np.diff(S)])
    df['cQ']=np.cumsum((Qin-Qout)*dt)
    return df

def GetHydrusData(mypath):
    RT=os.popen('tail -1 ../hydrus/%s/Balance.out'%mypath).read().split()[-1]
    print('%s runtime = ../hydrus/%s'%(mypath,RT))
    WB=pd.read_csv('../hydrus/%s/T_Level.out'%mypath,skiprows=[0,1,2,3,4,5,7],skipfooter=1,delim_whitespace=True,index_col=0,engine='python')
    WB.index.name=None
    Qin=-WB['vTop'].values*1000
    Qout=-WB['vBot'].values*1000
    S=WB['Volume'].values*1000
    start='1999-10-1'
    WB['t']=WB.index
    t=pd.to_datetime(start)+WB['t'].apply(pd.Timedelta, unit='D')
    dt=np.hstack([0,np.diff(WB.index)])
    
    df=pd.DataFrame(index=t)
    df['dt']=dt
    df['Qin']=Qin
    df['Qout']=Qout
    df['S']=S
    df['dS']=np.hstack([0,np.diff(S)])
    df['cQ']=np.cumsum((Qin-Qout)*dt)
    return df

def errplot(df,c='b',alpha=1.):
    err=df['dS']-df['cQ'].diff()
    pl.plot(df['Qin'],err,'.',color='%s'%c,alpha=alpha)

print('RUNTIME PERFORMANCE')
ire=GetIresonData('IresonRun') # ATS solution, rtol=1e-6
ire2=GetIresonData('IresonOptRun') # ATS solution, rtol=1e-7
hyd=GetHydrusData('InfiltrationProblem') # Hydrus 1D
celia=GetIresonData('celia') # Celia MPM solution

print('WATER BALANCE PERFORMANCE:')
print('IresonRun:  %.4f,  %.2e'%GetError(ire))
print('IresonRun2: %.4f,  %.2e'%GetError(ire2))
print('Celia:      %.4f,  %.2e'%GetError(celia))
print('Hydrus:     %.4f,  %.2e'%GetError(hyd))

#pl.figure()
#pl.plot(celia['S'],'.')
#pl.plot(ire['S'])
#pl.show()
#
#pl.figure()
#pl.plot(celia['Qout'],'.')
#pl.plot(ire['Qout'])
#pl.show()
#
#pl.figure()
#pl.plot(celia['Qin'],'.')
#pl.plot(ire['Qin'])
#pl.show()

pl.figure(figsize=(8,8))
pl.subplot(3,1,1)
pl.plot(hyd['dS']-hyd['cQ'].diff(),alpha=0.5,label='Hydrus 1D')
pl.plot(ire['dS']-ire['cQ'].diff(),label='ATS solution')
pl.plot(ire2['dS']-ire2['cQ'].diff(),alpha=0.75,label='ATS solution')
pl.plot(celia['dS']-celia['cQ'].diff(),alpha=0.75,label='celia model')
pl.gca().set_xticklabels('')
# pl.legend(fontsize=13)
pl.grid()
pl.ylabel('Timestep water \nbalance error (mm)',fontsize=13)
pl.subplot(3,1,2)
pl.plot(hyd['dS'].cumsum()-hyd['cQ'],alpha=0.5,label='Hydrus 1D')
pl.plot(ire['dS'].cumsum()-ire['cQ'],label='ATS solution, rtol=1E-6')
pl.plot(ire2['dS'].cumsum()-ire2['cQ'],alpha=0.75,label='ATS solution, rtol=1E-7')
pl.plot(celia['dS'].cumsum()-celia['cQ'],alpha=0.75,label='Celia MPM solution')
pl.legend(fontsize=13)
pl.grid()
pl.ylabel('Cumulated water \nbalance error (mm)',fontsize=13)

pl.subplot(3,2,5)
errplot(hyd,'tab:blue',0.5)
pl.xlabel('Infiltration (mm/d)')
pl.ylabel('Timestep water \nbalance error (mm)',fontsize=13)
pl.grid()

pl.subplot(3,2,6)
errplot(ire,'tab:orange',1.00)
errplot(ire2,'tab:green',0.75)
errplot(celia,'tab:red',0.5)
pl.xlabel('Infiltration (mm/d)')
pl.grid()

pl.subplots_adjust(left=0.18,right=0.96,top=0.96,bottom=0.1)
#pl.subplots_adjust(hspace=0.03,left=0.18,right=0.96,top=0.96,bottom=0.1)
pl.savefig('Figure06.png')

#pl.figure()
#errplot(ire,'r')
#errplot(ire2,'k')
#pl.legend(['Our model, rtol=1E-6','Our model, rtol=1E-7','Hydrus','SUMMA'])
#pl.grid()
#pl.show()
#
pl.figure()
pl.subplot(2,1,1)
pl.plot(hyd['S']['2004'],label='Hydrus')
pl.plot(ire['S']['2004'],label='ATS solution, rtol=1e-7')
pl.plot(ire2['S']['2004'],label='ATS solution, rtol=1e-6')
pl.plot(celia['S']['2004'],label='Celia',alpha=0.5)
pl.grid()
#pl.legend(ncol=3)
pl.gca().set_xticklabels([])
pl.ylabel('Storage (mm)')
pl.subplot(2,1,2)
pl.plot(hyd['Qout']['2004'],label='Hydrus')
pl.plot(ire['Qout']['2004'],label='ATS solution, rtol=1e-7')
pl.plot(ire2['Qout']['2004'],label='ATS solution, rtol=1e-6')
pl.plot(celia['Qout']['2004'],label='Celia',alpha=0.5)
pl.ylabel('Drainage (mm/d)')
pl.grid()
#pl.gca().set_yscale('log')
pl.subplots_adjust(hspace=0.04)
pl.legend()
pl.savefig('Figure07.png',dpi=300)
#pl.show()
#
