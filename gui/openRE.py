import dash
import dash_bootstrap_components as dbc
from dash import dcc, html
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np

# Import model libraries:
from Richards_FD import run_RE as run_RE_FD
from Richards_FP import run_RE as run_RE_FP
from Richards_ZF import run_RE as run_RE_ZF

from myFunctions import HydProps
from myFunctions import thetaFun
from myFunctions import CFun
from myFunctions import KFun
from myFunctions import runRE

# Options and defaults:
optionsLBC = [
    {'label': 'Free drainage', 'value': 'FD'},
    {'label': 'Zero flux', 'value': 'ZF'},
    {'label': 'Fixed psi', 'value': 'FP'}
]
defaultvalueLBC = 'FP'

optionsIC = [
    {'label': 'Fixed psi', 'value': 'FP'},
    {'label': 'Hydrostatic', 'value': 'HS'}
]
defaultvalueIC = 'HS'

# Dash app:
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])
app.title = "open RE"

input_style = {"width": "100px"}
label_style = {"marginRight": "10px"}


def make_input_row(label, id_, default):
    label=dcc.Markdown(label, mathjax=True)
    return dbc.Row([
        dbc.Col(dbc.Input(id=id_, type="text", value=default, style=input_style), width="auto"),
        dbc.Col(dbc.Label(label, html_for=id_, className="col-form-label"), width="auto")
    ], className="mb-2", align="center")


app.layout = dbc.Container([
    html.H1("Welcome to open RE - a simple Richard's Equation Solver", className="text-primary my-4"),

    html.H4("Soil hydraulic properties"),

    dbc.Row([
        dbc.Col([
            make_input_row("$\\psi_{min}$ (m):", "psi_min", "-10"),
            make_input_row("$\\psi_{max}$ (m):", "psi_max", "0."),
            make_input_row("$\\theta_R$ (-):", "thetaR", "0.1"),
            make_input_row("$\\theta_S$ (-):", "thetaS", "0.4"),
            make_input_row("$\\alpha$ (1/m):", "alpha", "0.5"),
            make_input_row("$n$:", "n", "1.8"),
            make_input_row("$K_S$ (m/d):", "KS", "0.2"),
            dbc.Button("Update Plot", id="update-button", color="primary", className="mt-2")
        ], width=3),

        dbc.Col([
            dcc.Graph(id='hyd-plot',mathjax=True)
        ], width=9)
    ], className="mb-5"),

    html.H4("Run openRE"),

    dbc.Row([
        dbc.Col([
            dbc.Label("Lower boundary condition:"),
            dcc.Dropdown(id='lowerBC', options=optionsLBC, value=defaultvalueLBC),

            dbc.Label("Initial condition:"),
            dcc.Dropdown(id='IC', options=optionsIC, value=defaultvalueIC),

            make_input_row("Times for infiltration pulses (d):", "tI", "0, 10,1000"),
            make_input_row("Infiltration pulses (m/d):", "Ipulses", "0.05, 0.05")
        ], width=6),

        dbc.Col([
            make_input_row("Runtime (d):", "RunTime", "10"),
            make_input_row("Time step (d):", "TimeStep", "0.25"),
            make_input_row("Soil depth (m):", "SoilDepth", "4"),
            make_input_row("Space step (m):", "SpaceStep", "0.1"),
            make_input_row("Lower/initial matric potential (m):", "psi_ini", "0."),
            dbc.Button("Run openRE", id="runRE-button", color="success", className="mt-2")
        ], width=6)
    ]),

    dcc.Graph(id='REplot', style={'display': 'none'}, mathjax=True),

    html.Footer("Andrew Ireson", style={'textAlign': 'center', 'padding': '10px', 'marginTop': '20px'})
], fluid=True)


@app.callback(
    Output('hyd-plot', 'figure'),
    Input('update-button', 'n_clicks'),
    State('psi_min', 'value'),
    State('psi_max', 'value'),
    State('thetaR', 'value'),
    State('thetaS', 'value'),
    State('alpha', 'value'),
    State('n', 'value'),
    State('KS', 'value')
)
def update_graph(n_clicks, psi_min, psi_max, thetaR, thetaS, alpha, n, KS):
    psi_min = float(psi_min)
    psi_max = float(psi_max)
    thetaR = float(thetaR)
    thetaS = float(thetaS)
    alpha = float(alpha)
    n = float(n)
    KS = float(KS)

    psi, theta, C, K = HydProps(psi_min, psi_max, thetaR, thetaS, alpha, n, KS)

    fig = make_subplots(
    rows=1, cols=3,
    subplot_titles=(r"$\theta(\psi)$", r"$C(\psi)$", r"$K(\psi)$")
    )

    fig.add_trace(go.Scatter(x=psi, y=theta, mode='lines', line=dict(color='blue')), row=1, col=1)
    fig.add_trace(go.Scatter(x=psi, y=C, mode='lines', line=dict(color='red')), row=1, col=2)
    fig.add_trace(go.Scatter(x=psi, y=K, mode='lines', line=dict(color='green')), row=1, col=3)

    fig.update_layout(showlegend=False, margin=dict(l=10, r=10, t=20, b=10))
    fig.update_xaxes(title_text=r'$\psi$', row=1, col=1)
    fig.update_yaxes(title_text=r'$\theta$', row=1, col=1)
    fig.update_xaxes(title_text=r'$\psi$', row=1, col=2)
    fig.update_yaxes(title_text=r'$C$', row=1, col=2)
    fig.update_xaxes(title_text=r'$\psi$', row=1, col=3)
    fig.update_yaxes(title_text=r'$K$', row=1, col=3)

    return fig


@app.callback(
    Output('REplot', 'figure'),
    Output('REplot', 'style'),
    Input('runRE-button', 'n_clicks'),
    State('RunTime', 'value'),
    State('TimeStep', 'value'),
    State('SoilDepth', 'value'),
    State('SpaceStep', 'value'),
    State('tI', 'value'),
    State('Ipulses', 'value'),
    State('psi_ini', 'value'),
    State('lowerBC', 'value'),
    State('IC', 'value')
)
def update_graph2(n_clicks, RunTime, TimeStep, SoilDepth, SpaceStep, tI, Ipulses, psi_ini, lowerBC, IC):
    if n_clicks is None:
        return {}, {'display': 'none'}
    t, z, psi, WB = runRE(RunTime, TimeStep, SoilDepth, SpaceStep, tI, Ipulses, psi_ini, lowerBC, IC)
    dt = float(TimeStep)

    fig = make_subplots(rows=1, cols=2, column_widths=[0.5, 0.5], subplot_titles=("Depth profile", "Time series"))
    i=0
    fig.add_trace(go.Scatter(x=psi[i, :], y=z, mode='lines', line=dict(color='blue'),name='Matric potential'), row=1, col=1)
    fig.add_trace(go.Scatter(x=psi[i, :] + float(SoilDepth) - z, y=z, mode='lines', line=dict(color='red'),name='Hydraulic head'), row=1, col=1)
    for i in range(1,len(t)):
        fig.add_trace(go.Scatter(x=psi[i, :], y=z, mode='lines', line=dict(color='blue'), showlegend=False), row=1, col=1)
        fig.add_trace(go.Scatter(x=psi[i, :] + float(SoilDepth) - z, y=z, mode='lines', line=dict(color='red'), showlegend=False), row=1, col=1)

    fig.add_trace(go.Scatter(x=WB.index, y=WB['S'] - WB['S'].iloc[0], mode='lines',name='Cum. change in storage'), row=1, col=2)
    fig.add_trace(go.Scatter(x=WB.index, y=WB['QIN'].cumsum() * dt - WB['QOUT'].cumsum() * dt, mode='markers',name='Cum. net flux'), row=1, col=2)
    fig.add_trace(go.Scatter(x=WB.index[1:], y=WB['QIN'].iloc[1:], mode='lines',name='Infiltration'), row=1, col=2)
    fig.add_trace(go.Scatter(x=WB.index[1:], y=WB['QOUT'].iloc[1:], mode='lines',name='Drainage'), row=1, col=2)

    fig.update_layout(xaxis_title='X', yaxis_title='Y', autosize=True,
                      margin=dict(l=10, r=10, t=10, b=10), 
                      yaxis=dict(range=[SoilDepth, 0]),
                      legend=dict(x=0.4,y=0.01,bordercolor="Black",  borderwidth=1)
)

    fig.update_xaxes(title_text='psi (m)', row=1, col=1)
    fig.update_yaxes(title_text='depth (m)', row=1, col=1)
    fig.update_xaxes(title_text='time (d)', row=1, col=2)
    fig.update_yaxes(title_text='Water balance', row=1, col=2)

    return fig, {'display': 'block'}

if __name__ == '__main__':
    app.run(debug=True, host="0.0.0.0")
