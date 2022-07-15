import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

def DoThePlot(runlist,dtlist,suflist,labels,fout,legpos,explode=False):
    cols=['m.','g.','bo','r.','c.'][:len(runlist)]
    pl.figure(figsize=(7,7))
    rc=0
    for run,l,suf,c in zip(runlist,labels,suflist,cols):
        err=np.zeros(len(dtlist))
        inQ=np.zeros(len(dtlist))
        for i,dt in enumerate(dtlist):
            dtstr=str(dt).replace('.','-')
            fname='%s/output/mb_dt%s.csv'%(run,dtstr)
            WB=pd.read_csv(fname,index_col=0)
            S=WB['S'].values
            if suf=='':
                Qin=WB['QIN'].values
                Qout=WB['QOUT'].values
            else:
                Qin=WB['QIN_%s'%suf].values
                Qout=WB['QOUT_%s'%suf].values
            dS=S[-1]-S[0]
            dQ=np.sum(Qin-Qout)*dt
            err[i]=1-dQ/dS
            inQ[i]=np.sum(Qin)*dt
        pl.subplot(2,1,1)
        pl.plot(dtlist,err,'-%s'%c,label=r'%s'%l)
#        if rc==2: 
#            yl1=pl.ylim()
#        elif rc>2:
#            pl.ylim(yl1)
        pl.subplot(2,1,2)
        pl.plot(dtlist,inQ,'-%s'%c,label=r'%s'%l)
#        if rc==2: 
#            yl2=pl.ylim()
#        elif rc>2:
#            pl.ylim(yl2)
        rc+=1

    pl.subplot(2,1,1)
    pl.legend(bbox_to_anchor=legpos,ncol=2,fontsize=13)
    pl.gca().set_xticklabels([])

    pl.grid()
    pl.ylabel('Global mass balance error',fontsize=13)
    #if explode:
    #    pl.ylabel(r'Global mass balance error $\times10^{-6}$')
    #else:
    #    pl.ylabel('Global mass balance error')
    #pl.ylim(0.75,1.25)

    pl.subplot(2,1,2)
    pl.xlabel('Time step (s)',fontsize=13)
    pl.ylabel('Cumulative inflow (cm)',fontsize=13)
    pl.grid()
    if fout=='MB_Celia': pl.ylim(2.2,2.5)
    pl.ylim(2.20,2.50)
    pl.subplots_adjust(left=0.13,hspace=0.03,right=0.97,top=0.82,bottom=0.1)
    pl.savefig('%s.png'%fout,dpi=300)


dtlist=np.array([0.1,1,3,10,20,30,40,60,90,120])

### DO CELIA MB PLOT
runlist=['../celiasolution/BackwardDiff','../celiasolution/Picard','../celiasolution/ModifiedPicard','../iresonsolution','../iresonsolution']
labels=['No iteration','Picard iteration','Celia MPM','ATS calculation step','ATS reporting step']
suflist=['BD','BD','BD','','BD']
fout='Figure03'
legpos=(0.05,1.03,0.95,0.1)

DoThePlot(runlist,dtlist,suflist,labels,fout,legpos)

