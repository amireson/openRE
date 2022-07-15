import numpy as np
import matplotlib.pyplot as pl

def DoThePlot(run_name,title,fout):
    dt=np.array([0.1, 1, 10, 30, 120])
    styles=['.','-','-','-','-']
    for dti,style in zip(dt,styles):
        dtstr=str(dti).replace('.','-')
        fname='%s/output/psi_dt%s.csv'%(run_name,dtstr)
        z,psi=np.loadtxt(fname,delimiter=',',skiprows=1,unpack=True)
        if dti<1:
            pl.plot(z,psi,style,label='dt=%.1f s'%dti)
        else:
            pl.plot(z,psi,style,label='dt=%.0f s'%dti)

    pl.ylim(-70,-10)
    pl.xlim(0,40)
    pl.grid()
    pl.xlabel('Depth (cm)',fontsize=13)
    pl.title(r'%s'%title,fontsize=13)

pl.figure(figsize=(12,4))
run_name='../celiasolution/BackwardDiff'
title='No iteration scheme'
fout='psi_BD'
pl.subplot(1,4,1)
DoThePlot(run_name,title,fout)
pl.ylabel('Pressure head (cm)',fontsize=13)

run_name='../celiasolution/Picard'
title='Picard iteration'
fout='psi_PI'
pl.subplot(1,4,2)
DoThePlot(run_name,title,fout)
pl.gca().set_yticklabels([])
xt=pl.gca().get_xticks()
pl.gca().set_xticks(xt[1:])

run_name='../celiasolution/ModifiedPicard'
title="Modified Picard Method"
fout='psi_MPM'
pl.subplot(1,4,3)
DoThePlot(run_name,title,fout)
xt=pl.gca().get_xticks()
pl.gca().set_xticks(xt[1:])
pl.gca().set_yticklabels([])

run_name='../iresonsolution'
title='ATS method'
fout='psi_MPM'
pl.subplot(1,4,4)
DoThePlot(run_name,title,fout)
pl.legend(fontsize=13)
xt=pl.gca().get_xticks()
pl.gca().set_xticks(xt[1:])
pl.gca().set_yticklabels([])

pl.subplots_adjust(wspace=0.03,left=0.06,right=0.97,top=0.92,bottom=0.15)
fout='Figure02'
pl.savefig('%s.png'%fout)

