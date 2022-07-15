import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

def GetError(df):
    dt=df.index[1]-df.index[0]
    Qin=df['QIN'].values*10
    Qout=df['QOUT'].values*10
    S=df['S'].values*10
    
    dS=S[-1]-S[0]
    dQ=np.sum((Qin-Qout)*dt)
    
    QIN=np.sum(Qin)*dt
    QOUT=np.sum(Qout)*dt
    err=1-dQ/dS
    
    MBerr=(dS-dQ)
    RMSE=np.sqrt(np.mean((np.diff(S)-(Qin[1:]-Qout[1:])*dt)**2))
    MBabserr=np.sum(np.abs(np.diff(S)-(Qin[1:]-Qout[1:])*dt))
    return MBerr, RMSE, QIN

names=['sand','siltloam','clay']
titles=['Hygiene sandstone','Silt Loam G.E. 3','Beit Netofa Clay']
dt=np.array([0.1,0.01,0.01])*60
pl.figure(figsize=(10,8))
for sb,fname in enumerate(['celiasolution','iresonsolution']):
    i=0
    for n in names:
        # Load mathias solution
        th,x1,x2,x3,x4,x5=np.loadtxt('../mathiassolution/%s.csv'%n,delimiter=',',unpack=True)
        # Load numerical solution
        x,th1,th2,th3,th4,th5=np.loadtxt('../%s/%s/theta.csv'%(fname,n),delimiter=',',unpack=True)
        # Load runtime:
        runtime=np.loadtxt('../%s/%s/runtime.csv'%(fname,n))
        WB=pd.read_csv('../%s/%s/mb.csv'%(fname,n),index_col=0,delimiter=',')
        MBerr,RMSE,QIN=GetError(WB)
        i+=1
        pl.subplot(2,3,i+3*sb)
        if sb==0: pl.title(titles[i-1],fontsize=14)
        #pl.plot(x1,th,'-',color='tab:blue',label='5 min')
        pl.plot(x,th1,'.',color='tab:blue',label='5 min')
        pl.plot(x1,th,'-k')
        #pl.plot(x2,th,'-',color='tab:orange',label='10 min')
        pl.plot(x,th2,'.',color='tab:orange',label='10 min')
        pl.plot(x2,th,'-k')
        #pl.plot(x3,th,'-',color='tab:green',label='20 min')
        pl.plot(x,th3,'.',color='tab:green',label='20 min')
        pl.plot(x3,th,'-k')
        #pl.plot(x4,th,'-',color='tab:red',label='50 min')
        pl.plot(x,th4,'.',color='tab:red',label='50 min')
        pl.plot(x4,th,'-k')
        #pl.plot(x5,th,'-',color='tab:purple',label='100 min')
        pl.plot(x,th5,'.',color='tab:purple',label='100 min')
        pl.plot(x5,th,'-k')
        if i==1: pl.plot(np.nan,np.nan,'-k',label='analytical solution')
        if i==1: pl.plot(np.nan,np.nan,'.k',label='numerical solution')
        if i==1 and sb==0: pl.ylabel("Celia's method\nWater content (-)",fontsize=13)
        if i==1 and sb==1: pl.ylabel("ATS method\nWater content (-)",fontsize=13)
        print(fname,n)
        print('dt = %.1f s'%dt[i-1])
        print('runtime = %.1f s'%runtime)
        print('QIN = %.1f mm'%QIN)
        print(('RMSE = %.1e mm'%RMSE).replace('e','E'))
        print(('MBerr = %.1e mm'%MBerr).replace('e','E'))
        #pl.text(2,0.12,'dt = %.1f s'%dt[i-1],fontsize=13)
        #pl.text(2,0.09,'runtime = %.1f s'%runtime,fontsize=13)
        #pl.text(2,0.06,'QIN = %.3f mm'%QIN,fontsize=13)
        #pl.text(2,0.03,'RMSE = %.3e mm'%RMSE,fontsize=13)
        pl.ylim(0,0.5)
        pl.xlim(0.1,100)
        pl.xscale('log')
        if sb==1: pl.xlabel('Distance (m)',fontsize=13)
        if sb==0: pl.gca().set_xticklabels([])
        if i==1 and sb==0: pl.legend(loc=2)
        if i!=1: pl.gca().set_yticklabels([])
        pl.grid()
        pl.subplots_adjust(left=0.10,bottom=0.1,right=0.97,top=0.95,wspace=0.06,hspace=0.06)
        xtl=pl.gca().get_xticks()
        #print(xtl)
        xtl=xtl[2:-2]
        if i==2: pl.gca().set_xticks(xtl)

pl.savefig('Figure05.png',dpi=300)
