import dash
from dash import html, dcc, Input, Output, State, clientside_callback
import numpy as np
import matplotlib.pyplot as plt
import io
import base64
import socket

class FastIonDiagnosticsApp:
    def __init__(self):
        # === Welcome message ===
        self.welcome_message = (
            "Welcome to the program!\n\n"
            "Explanation of the naming convention:\n"
            "For '21AEM.S8':\n"
            "  - 21:  \t   Module 2, Submodule 1\n"
            "  - AEM: \t   Type of the port\n"
            "  - S:   \t   NBI (FIDA) \n"
            "  - C:   \t   Gyrotron (CTS) \n"
            "  - 8:   \t   NBI or C number\n\n"
            "Matrix Size setting:\n"
            "This adjusts the size of each matrix.\n\n"
        )

        # === Initialize Dash ===
        self.app = dash.Dash(__name__)

        # === Build Layout ===
        self.app.layout = self.build_layout()

        # === Register Callbacks ===
        self.register_callbacks()

    @staticmethod
    def draw_matrix_plot(angle, matrix_size, delta_s, delta_j0, b0):
        """Generate matrix plot and return base64 image"""
        rows, cols = 7, 7
        cell_size = (matrix_size, matrix_size)

        fig, axs = plt.subplots(rows, cols, figsize=(7, 7))
        for i in range(rows):
            for j in range(cols):
                ax = axs[i, j]
                matrix = np.random.rand(*cell_size)
                ax.imshow(matrix, cmap='viridis', vmin=0, vmax=1)
                for spine in ax.spines.values():
                    spine.set_edgecolor('black')
                    spine.set_linewidth(1.0)
                ax.set_xticks([]), ax.set_yticks([]), ax.set_aspect('equal')

        plt.subplots_adjust(wspace=0, hspace=0)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        buf.seek(0)
        img_base64 = base64.b64encode(buf.read()).decode('utf-8')
        buf.close(), plt.close(fig)
        return f"data:image/png;base64,{img_base64}"

    @staticmethod
    def get_local_ip():
        """Get local network IP address for LAN link"""
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        try:
            s.connect(("8.8.8.8", 80))
            ip = s.getsockname()[0]
        except Exception:
            ip = "127.0.0.1"
        finally:
            s.close()
        return ip

    def build_layout(self):
        """Return the entire Dash layout"""
        return html.Div([
            html.Span("☰", id="open-sidebar-btn", className="open-btn"),

            # === Sidebar ===
            html.Div(
                id="sidebar",
                className="sidebar",
                children=[
                    html.Div([
                        html.Div([
                            html.H2([
                                html.Span("⚙️", className="settings-icon"),
                                " Settings"
                            ], className="sidebar-title-inline"),
                            html.Button("☰", id="close-sidebar-btn", className="close-btn")
                        ], className="sidebar-header"),

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
                                marks={10: "10", 30: "30", 50: "50", 70: "70", 100: "100"}
                            ),
                        ], className="slider-block"),

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
                        ], className="inputs-block"),

                        html.Button("Change Theme", id="toggle-theme-btn", n_clicks=0, className="toggle-btn"),

                    ], className="sidebar-content")
                ]
            ),

            html.H1("Fast Ion Diagnostics Placement Optimizer", className="main-title"),
            html.P("work in progress...", className="wip-label"),

            html.Div([
                html.Div([
                    html.P("Click the button to generate a random matrix grid."),
                    html.Button("Build Graph", id='build-graph-btn', n_clicks=0, className="build-graph-btn")
                ], className="settings-panel"),

                html.Div([
                    html.Img(id='matrix-graph', className='matrix-graph',
                             style={'max-width': '100%', 'height': 'auto', 'border-radius': '15px'})
                ], className="graph-panel"),

                html.Div([
                    html.Pre(
                        id="console-output",
                        children=self.welcome_message,
                        className="console-panel"
                    )
                ], className="console-container"),
            ], id="page-content", className="page-container"),

            dcc.Store(id='theme-store', data={'theme': 'light'}),
            html.Div(id='dummy-output-theme'),
            html.Div(id='dummy-output-sidebar')
        ])

    def register_callbacks(self):
        """Register all callbacks"""
        app = self.app

        # Theme switch + sidebar open/close
        clientside_callback(
        """
        function(data) {
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

            return '';
        }
        """,
        Output('dummy-output-theme', 'children'),
        Input('theme-store', 'data')
    )

    # === Clientside для выезда сайдбара ===
        clientside_callback(
        """
        function(openClicks, closeClicks) {
            const sidebar = document.getElementById('sidebar');
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
            return '';
        }
        """,
        Output('dummy-output-sidebar', 'children'),
        Input('open-sidebar-btn', 'n_clicks'),
        Input('close-sidebar-btn', 'n_clicks')
    )

        @app.callback(
            Output('theme-store', 'data'),
            Input('toggle-theme-btn', 'n_clicks'),
            State('theme-store', 'data'),
            prevent_initial_call=True
        )
        def toggle_theme(n_clicks, data):
            theme = data['theme']
            return {'theme': 'dark' if theme == 'light' else 'light'}

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
            return self.draw_matrix_plot(angle, matrix_size, delta_s, delta_j0, b0)

    def run(self):
        """Run server with local IP info"""
        host = "0.0.0.0"
        port = 8050
        local_ip = self.get_local_ip()
        print("=" * 50)
        print(f" 🚀 Your Dash app is running at:")
        print(f" 👉 Local:   http://127.0.0.1:{port}")
        print(f" 👉 Network: http://{local_ip}:{port}")
        print("=" * 50)
        self.app.run(host=host, port=port, debug=True)


if __name__ == '__main__':
    app = FastIonDiagnosticsApp()
    app.run()
