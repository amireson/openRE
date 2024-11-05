# %load openRE.py
# %load openRE.py
import dash
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

from Richards_FD import run_RE as run_RE_FD
from Richards_FP import run_RE as run_RE_FP
from Richards_ZF import run_RE as run_RE_ZF

def HydProps(psi_min,psi_max,thetaR, thetaS,alpha,n,KS):
    global pars
    psi=np.linspace(psi_min,psi_max)
    pars={}
    pars['thetaR']=thetaR
    pars['thetaS']=thetaS
    pars['alpha']=alpha
    pars['n']=n
    pars['m']=1-1/n
    pars['Ks']=KS
    pars['neta']=0.5
    pars['Ss']=1e-6

    theta=thetaFun(psi,pars)
    C=CFun(psi,pars)
    K=KFun(psi,pars)

    return psi,theta,C,K

def thetaFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se

def CFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
    return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh

def KFun(psi,pars):
    Se=(1+(psi*-pars['alpha'])**pars['n'])**(-pars['m'])
    Se[psi>0.]=1.0
    return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2

def runRE(RunTime,TimeStep,SoilDepth,SpaceStep,tI,Ipulses,psi_ini,lowerBC,IC):
    
    global pars
    
    # Time grid:
    tN=float(RunTime)

    # Spatial grid:
    dz=float(SpaceStep)
    zN=float(SoilDepth)
    z=np.arange(dz/2,zN,dz)
    n=len(z)
    # z=np.hstack([0,z,zN])
    #z=z[-1]-z

    # Initial condition:
    if IC=='HS':
        psi0=z-zN
    elif IC=='FP':
        psi0=np.zeros(n)+float(psi_ini)
    else:
        print('error')

    
    
    psiB=float(psi_ini)
    
    dt=float(TimeStep)
    t=np.arange(0,tN+dt,dt)
    nt=len(t)

    # Boundary conditions:
    tI=[float(i) for i in tI.split(',')]
    Ipulses=[float(i) for i in Ipulses.split(',')]

    I=np.zeros(len(t))
    c=0

    for ti in range(len(t)):
        if t[ti]>=tI[c+1]:
            c+=1
        I[ti]=Ipulses[c]
    
    
    # I=(0.5-np.cos(t*2*np.pi/365)/2)*(maxInf-minInf)+minInf
    
    BC_T=I+np.zeros(nt)
    BC_B=psiB+np.zeros(nt)
    
    if lowerBC=='FD':
        psi,WB,runtime=run_RE_FD(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)
    elif lowerBC=='FP':
        psi,WB,runtime=run_RE_FP(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)
    elif lowerBC=='ZF':    
        psi,WB,runtime=run_RE_ZF(dt,t,dz,zN,n,psi0,BC_T,BC_B,pars)
        
    return t,z,psi

optionsLBC=[{'label': 'Free drainage', 'value': 'FD'},
         {'label': 'Zero flux', 'value': 'ZF'},
         {'label': 'Fixed psi', 'value': 'FP'}]
defaultvalueLBC='FD'

optionsIC=[{'label': 'Fixed psi', 'value': 'FP'},
         {'label': 'Hydrostatic', 'value': 'HS'}]
defaultvalueIC='FP'

app = dash.Dash(__name__)

app.title = "open RE"

app.layout = html.Div([

    html.Header("Welcome to open RE - a simple Richard's Equation Solver", style={'textAlign': 'center', 'fontSize': '26px', 'padding': '10px'}),
    
    html.Div([
        html.P("Soil hydraulic properties"),
    ],style={'textAlign': 'left', 'fontSize': '22px', 'padding': '10px'}),
    
    html.Div([
        
    html.Div([
        
    # Text boxes for soil hydraulic properties
    dcc.Input(
        id='psi_min', 
        type='text', 
        value='-10',
        style={'width': '75px'}
    ),
    html.Label(" psi_min (m): "),
    html.Br(),
    dcc.Input(
        id='psi_max', 
        type='text', 
        value='0.',
        style={'width': '75px'}
    ),
    html.Label(" psi_max (m): "),

    html.Br(),

    dcc.Input(
        id='thetaR', 
        type='text', 
        value='0.1',
        style={'width': '75px'}
    ),
    html.Label(" theta_R: "),
    html.Br(),
    dcc.Input(
        id='thetaS', 
        type='text', 
        value='0.4', 
        style={'width': '75px'}
    ),
    html.Label(" theta_S: "),
    html.Br(),
    dcc.Input(
        id='alpha', 
        type='text', 
        value='0.5', 
        style={'width': '75px'}
    ),
    html.Label(" alpha (1/m): "),
    html.Br(),
    dcc.Input(
        id='n', 
        type='text', 
        value='1.8', 
        style={'width': '75px'}
    ),
    html.Label(" n: "),
    html.Br(),
    dcc.Input(
        id='KS', 
        type='text', 
        value='0.2', 
        style={'width': '75px'}
    ),
    html.Label(" K_S (m/d): "),
    html.Br(),
    html.Br(),
    html.Button('Update Plot', id='update-button', n_clicks=0),
    html.Br(),
    ], style={'width': '18%', 'display': 'inline-block', 'padding': '10px'}),

    html.Div([
    # Plot hydraulic properties
    
    dcc.Graph(id='hyd-plot'),

    ], style={'width': '78%', 'display': 'inline-block', 'padding': '10px'}),
    ], style={'display': 'flex', 'justifyContent': 'center'}),

    html.Div([
        html.P("Run openRE"),
    ],style={'textAlign': 'left', 'fontSize': '22px', 'padding': '10px'}),
    
    html.Div([

    html.Div([
    # Options for model run  
    html.Label(" Lower boundary condition: "),
    dcc.Dropdown(
        id='lowerBC',
        options=optionsLBC,
        value=defaultvalueLBC,
        style={'width': '300px'}),

    html.Label(" Initial condition: "),
    dcc.Dropdown(
        id='IC',
        options=optionsIC,
        value=defaultvalueIC,
        style={'width': '300px'}),
    html.Br(),
    html.Br(),

    dcc.Input(
        id='tI', 
        type='text', 
        value='0, 10,1000', 
        style={'width': '150px'}
    ),
    html.Label(" Times for infiltration pulses (d): "),

    html.Br(),
    dcc.Input(
        id='Ipulses', 
        type='text', 
        value='0.1, 0.1', 
        style={'width': '150px'}
    ),
    html.Label(" Infiltration pulses (m/d): "),
    
    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px'}),

    html.Div([    
    # Text boxes for model run
    html.Br(),
    dcc.Input(
        id='RunTime', 
        type='text', 
        value='10', 
        style={'width': '75px'}
    ),
    html.Label(" Runtime (d): "),

    html.Br(),
    dcc.Input(
        id='TimeStep', 
        type='text', 
        value='0.25', 
        style={'width': '75px'}
    ),
    html.Label(" Time step (d): "),
    
    html.Br(),
    dcc.Input(
        id='SoilDepth', 
        type='text', 
        value='4', 
        style={'width': '75px'}
    ),
    html.Label(" Soil depth (m): "),

    html.Br(),
    dcc.Input(
        id='SpaceStep', 
        type='text', 
        value='0.1', 
        style={'width': '75px'}
    ),
    html.Label(" Space step (m): "),

    html.Br(),
    dcc.Input(
        id='psi_ini', 
        type='text', 
        value='-5', 
        style={'width': '75px'}
    ),
    html.Label(" Lower/initial matric potential (m): "),

    html.Br(),
    html.Br(),
    html.Button('Run openRE', id='runRE-button', n_clicks=0),
    ], style={'width': '48%', 'display': 'inline-block', 'padding': '10px'}),

    html.Div([    
    # Run model and plot output
    dcc.Graph(id='REplot', style={'display': 'none'})
    ]),
    ]),
    
    html.Footer("Andrew Ireson", style={'textAlign': 'right', 'padding': '10px', 'marginTop': '20px'})
    
], style={'padding': '20px'})

@app.callback(
    Output('hyd-plot', 'figure'),
    [Input('update-button', 'n_clicks')],
    [State('psi_min', 'value'),
     State('psi_max', 'value'),
     State('thetaR', 'value'),
     State('thetaS', 'value'),
     State('alpha', 'value'),
     State('n', 'value'),
     State('KS', 'value')]
)
def update_graph(n_clicks, psi_min,psi_max,thetaR, thetaS,alpha,n,KS):
    psi_min = float(psi_min)
    psi_max = float(psi_max)
    thetaR = float(thetaR)
    thetaS = float(thetaS)
    alpha = float(alpha)
    n = float(n)
    KS = float(KS)

    psi,theta,C,K=HydProps(psi_min,psi_max,thetaR, thetaS,alpha,n,KS)

    # Create a figure with 3 subplots
    fig = make_subplots(rows=1, cols=3, subplot_titles=("theta-psi", "C-psi", "K-psi"))

    # Add each plot to the figure
    fig.add_trace(go.Scatter(x=psi, y=theta, mode='lines', name='sin(x)', line=dict(color='blue')), row=1, col=1)
    fig.add_trace(go.Scatter(x=psi, y=C, mode='lines', name='cos(x)', line=dict(color='red')), row=1, col=2)
    fig.add_trace(go.Scatter(x=psi, y=K, mode='lines', name='tan(x)', line=dict(color='green')), row=1, col=3)

    # Update layout
    fig.update_layout(showlegend=False,margin=dict(l=0, r=0, t=0, b=0))

    # Axes labels
    fig.update_xaxes(title_text='psi', row=1, col=1)
    fig.update_yaxes(title_text='theta', row=1, col=1)
    fig.update_xaxes(title_text='psi', row=1, col=2)
    fig.update_yaxes(title_text='C', row=1, col=2)
    fig.update_xaxes(title_text='psi', row=1, col=3)
    fig.update_yaxes(title_text='K', row=1, col=3)

    return fig

@app.callback(
    Output('REplot', 'figure'),
    Output('REplot', 'style'),  # Control the visibility
    Input('runRE-button', 'n_clicks'),
    [State('RunTime', 'value'),
     State('TimeStep', 'value'),
     State('SoilDepth', 'value'),
     State('SpaceStep', 'value'),
     State('tI', 'value'),
     State('Ipulses', 'value'),
     State('psi_ini', 'value'),
     State('lowerBC', 'value'),
     State('IC', 'value')]
)
def update_graph2(n_clicks,RunTime,TimeStep,SoilDepth,SpaceStep,tI,Ipulses,psi_ini,lowerBC,IC):
    global pars
    if n_clicks > 0:
        # Generate random data for the second graph
        t,z,psi=runRE(RunTime,TimeStep,SoilDepth,SpaceStep,tI,Ipulses,psi_ini,lowerBC,IC)

        # Create the figure for the second graph
        fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], subplot_titles=("Depth profile", "Time series"))
        # fig = go.Figure()
        for i in range(len(t)):
            fig.add_trace(go.Scatter(x=psi[i,:], y=z, mode='lines', line=dict(color='blue')),row=1,col=1)
            fig.add_trace(go.Scatter(x=psi[i,:]+float(SoilDepth)-z, y=z, mode='lines', line=dict(color='red')),row=1,col=1)

        for i in range(len(z)):
            fig.add_trace(go.Scatter(x=t, y=psi[:,i], mode='lines'),row=1,col=2)
            
        fig.update_layout(xaxis_title='X', yaxis_title='Y', 
                            width=None,
                            autosize=True,
                            margin=dict(l=10, r=10, t=10, b=10),
                            showlegend=False, yaxis=dict(range=[SoilDepth, 0]))

         # Axes labels
        fig.update_xaxes(title_text='psi (m)', row=1, col=1)
        fig.update_yaxes(title_text='depth (m)', row=1, col=1)
        fig.update_xaxes(title_text='time (d)', row=1, col=2)
        fig.update_yaxes(title_text='psi (m)', row=1, col=2)
        
        
        return fig, {'display': 'block'}  # Show the graph
    else:
        return {}, {'display': 'none'}  # Keep it hidden

# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
