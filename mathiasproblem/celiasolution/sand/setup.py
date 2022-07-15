pars={}
pars['thetaR']=0.153
pars['thetaS']=0.250
pars['alpha']=0.0079
pars['n']=10.4 
pars['m']=1-1/pars['n']
pars['Ks']=108.
pars['neta']=0.5
pars['Ss']=0.

# Time grid:
dt=0.1/24./60.
tN=100./24./60.

# Spatial grid:
dz=0.25
zN=100.

# Boundary conditions:
Se0=0.01
SeI=0.99
psiT=-(SeI**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
psiB=-(Se0**(-1/pars['m'])-1)**(1/pars['n'])/pars['alpha']
