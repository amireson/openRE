import pandas as pd
import numpy as np
import matplotlib.pyplot as pl

qT=pd.read_csv('infiltration.dat',index_col=0,parse_dates=True)
P=qT.iloc[:,0].values
print(qT)

#for y in range(1999,2010):
#    s='%d-10'%y
#    e='%d-09'%(y+1)
#    pl.plot(qT[s:e].iloc[:,0].cumsum())
#pl.show()

Pann=qT.iloc[:,0].resample('A-SEP').sum().values

print('Max P: %.0f'%np.max(P))
print('Mean annual P: %.0f'%np.mean(Pann))
print('Min annual P: %.0f'%np.min(Pann))
print('Max annual P: %.0f'%np.max(Pann))

print(qT[P>40])
