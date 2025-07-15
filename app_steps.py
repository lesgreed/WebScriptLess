import dash
from dash import html, dcc, Input, Output, State, clientside_callback
import numpy as np
import matplotlib.pyplot as plt
import io
import base64

# ==========================================================
# 🔵 UTILS: Function to generate a 7x7 grid of random matrices
# Each cell is 10x10, borders added for clear separation.
# ==========================================================

def draw_matrix_plot(angle, matrix_size, delta_s, delta_j0, b0):
    rows, cols = 7, 7
    cell_size = (matrix_size, matrix_size)

    fig, axs = plt.subplots(rows, cols, figsize=(7, 7))

    for i in range(rows):
        for j in range(cols):
            ax = axs[i, j]
            matrix = np.random.rand(*cell_size)
            ax.imshow(matrix, cmap='viridis', vmin=0, vmax=1)

            # Add border around each matrix
            for spine in ax.spines.values():
                spine.set_edgecolor('black')
                spine.set_linewidth(1.0)

            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_aspect('equal')

    plt.subplots_adjust(wspace=0, hspace=0)  # No spacing between cells

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)
    return f"data:image/png;base64,{img_base64}"


# ==========================================================
# 🟢 INITIALIZE DASH APP
# ==========================================================

app = dash.Dash(__name__)

# ==========================================================
# 🟢 APP LAYOUT: Sidebar, title, config panel, and output canvas
# ==========================================================

app.layout = html.Div([

    # ☰ button to open sidebar
    html.Span("☰", id="open-sidebar-btn", className="open-btn"),

    # === Sidebar ===
html.Div(
    id="sidebar",
    className="sidebar",
    children=[
        html.Div([

            # Sidebar header
            html.Div([
                html.H2([
                    html.Span("⚙️", className="settings-icon"),
                    " Settings"
                ], className="sidebar-title-inline"),
                html.Button("☰", id="close-sidebar-btn", className="close-btn")
            ], className="sidebar-header"),

 # === Sliders go immediately under header ===
    html.Div([
        html.Label("Angle", className="label"),
        dcc.Slider(
            id="angle-slider",
            min=0, max=90, step=1, value=90,
            tooltip={"always_visible": False, "placement": "bottom"},
            marks={0: "0°", 90: "90°"}
        ),

        html.Label("Matrix Size", className="label"),
        dcc.Slider(
            id="matrix-size-slider",
            min=10, max=100, step=5, value=10,
            tooltip={"always_visible": False, "placement": "bottom"},
            marks={10: "10",30: "30", 50: "50", 70: "70", 100: "100"}
        ),
    ], className="slider-block"),


# === Inputs Block Below Sliders ===
html.Div([

    html.Div([

    html.Div([
        html.Label("Δs =", className="inline-label"),
        dcc.Input(
            id="delta-s-input",
            type="number",
            value=0.05,
            step=0.01,
            className="sidebar-input"
        ),
    ], className="input-pair"),

    html.Div([
        html.Label("ΔJ₀ =", className="inline-label"),
        dcc.Input(
            id="delta-j0-input",
            type="number",
            value=0.005,
            step=0.001,
            className="sidebar-input"
        ),
    ], className="input-pair"),

    html.Div([
        html.Label("B₀ =", className="inline-label"),
        dcc.Input(
            id="b0-input",
            type="number",
            value=2.52,
            step=0.01,
            className="sidebar-input"
        ),
    ], className="input-pair"),

], className="inputs-block")


], className="slider-block"),



    # === Theme toggle button at the bottom ===
    html.Button("Change Theme", id="toggle-theme-btn", n_clicks=0, className="toggle-btn"),

        ], className="sidebar-content")
    ]
),
    # === Main Title ===
    html.H1("Fast Ion Diagnostics Placement Optimizer", className="main-title"),

    html.P("work in progress...", className="wip-label"),

    # === Main Panel: Config + Graph ===
    html.Div([
        html.Div([
            html.H3("Configuration Panel"),
            html.P("Click the button to generate a random matrix grid."),
            html.Button("Build Graph", id='build-graph-btn', n_clicks=0, className="build-graph-btn")
        ], className="settings-panel"),
        

        html.Div([
            html.H3("Graph Output"),
            html.Img(id='matrix-graph', className='matrix-graph', style={'max-width': '100%','height': 'auto',  'border-radius': '15px'})
        ], className="graph-panel")
    ], id="page-content",className="page-container"),

    # === Theme store and dummy outputs ===
    dcc.Store(id='theme-store', data={'theme': 'light'}),
    html.Div(id='dummy-output-theme'),
    html.Div(id='dummy-output-sidebar')
])



html.Div([
    html.H3("Configuration Panel"),

    html.P("Click the button to generate a random matrix grid."),

    # === Delta s ===
    html.Label("Δs:"),
    dcc.Input(
        id="delta-s-input",
        type="number",
        value=0.05,  # Default value
        step=0.01,
        style={"width": "100%", "margin-bottom": "10px"}
    ),

    # === Delta J0 ===
    html.Label("ΔJ0:"),
    dcc.Input(
        id="delta-j0-input",
        type="number",
        value=0.005,
        step=0.001,
        style={"width": "100%", "margin-bottom": "10px"}
    ),

    # === B0 ===
    html.Label("B0:"),
    dcc.Input(
        id="b0-input",
        type="number",
        value=2.52,
        step=0.01,
        style={"width": "100%", "margin-bottom": "10px"}
    ),

    # === Your Build Button ===
    html.Button("Build Graph", id='build-graph-btn', n_clicks=0, className="build-graph-btn")
], className="settings-panel"),




# ==========================================================
# 🟢 Clientside callback: Apply theme (light or dark)
# ==========================================================

clientside_callback(
    """
    function(data, openClicks, closeClicks) {
        // Theme switcher
        const theme = data.theme || 'light';
        document.body.className = theme + "-theme";

        const sidebar = document.getElementById('sidebar');
        if (theme === 'light') {
            sidebar.classList.remove('dark-sidebar');
            sidebar.classList.add('light-sidebar');
        } else {
            sidebar.classList.remove('light-sidebar');
            sidebar.classList.add('dark-sidebar');
        }

        // Sidebar shift logic
        const pageContent = document.getElementById('page-content');
        if (openClicks === undefined) openClicks = 0;
        if (closeClicks === undefined) closeClicks = 0;

        if (openClicks > closeClicks) {
            sidebar.style.width = '250px';
            pageContent.style.marginLeft = '250px';
        } else {
            sidebar.style.width = '0';
            pageContent.style.marginLeft = '0';
        }

        return ['', ''];
    }
    """,
    Output('dummy-output-theme', 'children'),
    Input('theme-store', 'data'),
    Input('open-sidebar-btn', 'n_clicks'),
    Input('close-sidebar-btn', 'n_clicks')
)


# ==========================================================
# 🟢 Python callback: Toggle theme and update store
# ==========================================================

@app.callback(
    Output('theme-store', 'data'),
    Input('toggle-theme-btn', 'n_clicks'),
    State('theme-store', 'data'),
    prevent_initial_call=True
)
def toggle_theme(n_clicks, data):
    theme = data['theme']
    new_theme = 'dark' if theme == 'light' else 'light'
    return {'theme': new_theme}


# ==========================================================
# 🟢 Clientside callback: Open / Close Sidebar
# ==========================================================

clientside_callback(
    """
    function(openClicks, closeClicks) {
        const sidebar = document.getElementById('sidebar');
        if (openClicks === undefined) openClicks = 0;
        if (closeClicks === undefined) closeClicks = 0;
        if (openClicks > closeClicks) {
            sidebar.style.width = '250px';
        } else {
            sidebar.style.width = '0';
        }
        return '';
    }
    """,
    Output('dummy-output-sidebar', 'children'),
    Input('open-sidebar-btn', 'n_clicks'),
    Input('close-sidebar-btn', 'n_clicks')
)

# ==========================================================
# 🟢 Python callback: Build matrix graph when button clicked
# ==========================================================

@app.callback(
    Output('matrix-graph', 'src'),
    Input('build-graph-btn', 'n_clicks'),
    State('angle-slider', 'value'),
    State('matrix-size-slider', 'value'),
    State('delta-s-input', 'value'),
    State('delta-j0-input', 'value'),
    State('b0-input', 'value'),
    prevent_initial_call=True
)
def build_matrix_graph(n_clicks, angle, matrix_size, delta_s, delta_j0, b0):
    # 👇 Теперь ты передаёшь ВСЁ
    return draw_matrix_plot(angle, matrix_size, delta_s, delta_j0, b0)


# ==========================================================
# 🟢 Run the server
# ==========================================================

if __name__ == '__main__':
    app.run(debug=True)
