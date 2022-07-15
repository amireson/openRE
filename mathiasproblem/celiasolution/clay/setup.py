pars={}
pars['thetaR']=0.
pars['thetaS']=0.446
pars['alpha']=0.00152
pars['n']=1.17 
pars['m']=1-1/pars['n']
pars['Ks']=0.082
pars['neta']=0.5
pars['Ss']=0.

# Time grid:
dt=0.01/24./60.
tN=100./24./60.

# Spatial grid:
dz=0.0025
zN=1.

# Boundary conditions:
Se0=0.01
SeI=0.99
psiT=-(SeI**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
psiB=-(Se0**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
