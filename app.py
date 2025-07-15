# app.py
import dash
from dash import html, dcc, Input, Output, State
import plotly.express as px
import os
from datetime import datetime

from modules import Data
from modules import calculus


app = dash.Dash(__name__)
data_instance = Data()
Bget = calculus()


Name_Ports = []
Name_NBI = []
diagnostics = []

# Layout
app.layout = html.Div([
    html.H1("W7-X Checking (Dash version)"),

    html.Div([
        html.Label("Appearance Mode:"),
        dcc.Dropdown(
            id='appearance-mode',
            options=[
                {'label': 'Light', 'value': 'Light'},
                {'label': 'Dark', 'value': 'Dark'},
                {'label': 'System', 'value': 'System'}
            ],
            value='Dark'
        ),

        html.Br(),

        html.Label("Viewing Angle:"),
        dcc.Slider(30, 90, 1, value=90, id='angle-slider'),
        html.Div(id='angle-output'),

        html.Br(),

        html.Label("Matrix Size:"),
        dcc.Slider(1, 10, 1, value=1, id='matrix-slider'),
        html.Div(id='matrix-output'),

        html.Br(),

        html.Label("Δs:"),
        dcc.Input(id='delta-s-input', value=0.05, type='number'),

        html.Label("ΔJ0 %:"),
        dcc.Input(id='delta-j0-input', value=0.5, type='number'),

        html.Label("B0:"),
        dcc.Input(id='B0-input', value=2.520, type='number'),

        html.Br(),

        html.Label("Select NBI:"),
        dcc.Dropdown(id='nbi-dropdown',
                     options=[{'label': f'NBI_{i}', 'value': f'NBI_{i}'} for i in range(1, 9)],
                     value='NBI_1'),

        html.Label("Select Port:"),
        dcc.Dropdown(id='port-dropdown'),

        html.Label("Diagnostic type:"),
        dcc.Dropdown(
            id='diagnostic-dropdown',
            options=[{'label': 'FIDA', 'value': 'FIDA'}, {'label': 'CTS', 'value': 'CTS'}],
            value='FIDA'
        ),

        html.Button('Add Port', id='add-port-button'),
        html.Button('Build', id='build-button'),

        html.Div(id='log-box', style={'whiteSpace': 'pre-line', 'border': '1px solid #ccc',
                                      'padding': '10px', 'marginTop': '10px'})
    ], style={'width': '30%', 'display': 'inline-block', 'verticalAlign': 'top'}),

    html.Div([
        dcc.Graph(id='matrix-graph')
    ], style={'width': '65%', 'display': 'inline-block', 'padding': '20px'})
])

# Update sliders
@app.callback(Output('angle-output', 'children'), Input('angle-slider', 'value'))
def update_angle_output(value):
    return f'Angle: {value}°'

@app.callback(Output('matrix-output', 'children'), Input('matrix-slider', 'value'))
def update_matrix_output(value):
    return f'Matrix Size: {value*10}'

# Update port dropdown when NBI changes
@app.callback(
    Output('port-dropdown', 'options'),
    Input('nbi-dropdown', 'value'),
)
def update_port_options(selected_nbi):
    index = int(selected_nbi.split('_')[1]) - 1
    ports = data_instance.port_for_nbi(index, 90, 10)[3]
    return [{'label': p, 'value': p} for p in ports]

# Add port and log
@app.callback(
    Output('log-box', 'children'),
    Input('add-port-button', 'n_clicks'),
    State('nbi-dropdown', 'value'),
    State('port-dropdown', 'value'),
    State('diagnostic-dropdown', 'value'),
    State('log-box', 'children'),
    prevent_initial_call=True
)
def add_port(n_clicks, selected_nbi, selected_port, diagnostic_type, log):
    Name_NBI.append(selected_nbi)
    Name_Ports.append(selected_port)
    diagnostics.append(diagnostic_type)
    timestamp = datetime.now().strftime("%H:%M:%S")
    return f"{log}\n[{timestamp}] Added Port: {selected_port}, NBI: {selected_nbi}, Type: {diagnostic_type}"

# Build graph
@app.callback(
    Output('matrix-graph', 'figure'),
    Input('build-button', 'n_clicks'),
    State('angle-slider', 'value'),
    State('matrix-slider', 'value'),
    State('delta-s-input', 'value'),
    State('delta-j0-input', 'value'),
    State('B0-input', 'value'),
    prevent_initial_call=True
)
def build_graph(n_clicks, angle, scale, delta_s, delta_j0, B0):
    # Пример: здесь вызывай свои методы
    # Например:
    # all_results = data_instance.data_already_input(...)
    # Матрицы = Bget.compute_matrix(...)
    # Для примера строим фиктивную картинку:
    fig = px.imshow([[1, 2], [3, 4]], color_continuous_scale='jet')
    fig.update_layout(title=f'Matrix @ angle={angle}, scale={scale*10}')
    return fig

if __name__ == '__main__':
    app.run(debug=True)
