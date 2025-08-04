import dash
from dash import html, dcc, Input, Output, State, clientside_callback
import numpy as np
import matplotlib.pyplot as plt
import io
import base64
import socket
from modules import Data
from modules import calculus
import os
import customtkinter as ctk
import threading
import time
from dash import callback_context
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from flask_caching import Cache
from dash import ctx 
from datetime import datetime
import tempfile
import os



# ==========================================================
# 🟢 UTILS
# ==========================================================


#-----------------------directory---------------------------

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(BASE_DIR, "config_files")
calc_lock = threading.Lock()

# --------------------- server info -----------------------
def get_local_ip():
    """Get local IP for sharing the app on local network."""
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
    except Exception:
        ip = "127.0.0.1"
    finally:
        s.close()
    return ip

# -------------- Button config -----------------
def get_config_files():
    files = []
    for filename in os.listdir(CONFIG_DIR):
        if filename.endswith(".txt"):
            files.append(os.path.join(CONFIG_DIR, filename))
    return files


from flask import request

import uuid
from flask import g

_session_ids = {}  # локальное хранилище для пользователей (по IP)

def user_cache_key(base_key="user_result_array"):
    """Создаёт ключ для кэша, который сохраняется, пока пользователь не обновит страницу."""
    try:
        user_ip = request.remote_addr or "anon"
    except RuntimeError:
        user_ip = "server"

    # ✅ Если для этого пользователя ещё нет session_id → создаём
    if user_ip not in _session_ids:
        _session_ids[user_ip] = str(uuid.uuid4())

    session_id = _session_ids[user_ip]
    return f"{base_key}_{session_id}"


# ==========================================================
# 🟢 CONSTANTS
# ==========================================================

welcome_message = (
    "Welcome!\n\n"
    "Explanation of the naming convention:\n"
    "For '21AEM.S8':\n"
    "  - 21:  Module 2, Submodule 1\n"
    "  - AEM: Type of the port\n"
    "  - S:   NBI (FIDA)\n"
    "  - C:   Gyrotron (CTS)\n"
    "  - 8:   NBI or Gyrotron number\n\n"
    "Matrix Size setting:\n"
    "This adjusts the size of each matrix.\n"
    "-------------------------------------\n"
)

# ==========================================================
# 🟢 INITIALIZE 
# ==========================================================

app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,
    prevent_initial_callbacks=True
)
app.enable_duplicate_callbacks = True
cache = Cache(app.server, config={ 'CACHE_TYPE': 'filesystem', 'CACHE_DIR': 'cache-dir','CACHE_DEFAULT_TIMEOUT': 3600})
data_instance = Data()
calc = calculus()


# ==========================================================
# 🟢 APP LAYOUT
# ==========================================================
app.layout = html.Div([
    #-------------------------------------------------variables----------------------------------------------------
    dcc.Store(id="calculation-flag", data=False),
    dcc.Store(id="log-list-store", data=[[welcome_message]]),    
    dcc.Store(id="log-store", data=welcome_message),   
    dcc.Store(id="precalc-flag", data=False),
    dcc.Store(id="standard-flag", data=False),
    dcc.Store(id="add-port-flag", data=False),
    dcc.Store(id="save-results-flag", data=False),


     
    # ☰ Open sidebar button
    html.Span("☰", id="open-sidebar-btn", className="open-btn"),
    
    # === Sidebar ===
    html.Div( id="sidebar",className="sidebar", children=[
            html.Div([ 
                
            html.Div([
                    html.H2([html.Span("⚙️"), " Settings"], className="sidebar-title-inline"),
                    html.Button("☰", id="close-sidebar-btn", className="close-btn") 
                     ], className="sidebar-header"),

            html.Div([
                    html.Label("Angle", className="label"),
                    dcc.Slider(id="angle-slider", min=0, max=90, step=1, value=90,tooltip={"always_visible": False, "placement": "bottom"},marks={0: "0°", 90: "90°"}),
                    html.Label("Matrix Size", className="label"),
                    dcc.Slider(id="matrix-size-slider",min=10, max=100, step=5, value=10,tooltip={"always_visible": False, "placement": "bottom"},marks={10: "10", 50: "50", 100: "100"}),
                     ], className="slider-block"),

            html.Div([
                    html.Div([
                        html.Label("Δs =", className="inline-label"),
                        dcc.Input(id="delta-s-input", type="number", value=0.05, step=0.01, className="sidebar-input")
                            ], className="input-pair"),

                    html.Div([
                        html.Label("ΔJ₀ =", className="inline-label"),
                        dcc.Input(id="delta-j0-input", type="number", value=0.005, step=0.001, className="sidebar-input")
                             ], className="input-pair"),

                    html.Div([
                        html.Label("B₀ =", className="inline-label"),
                        dcc.Input(id="b0-input", type="number", value=2.52, step=0.01, className="sidebar-input")
                             ], className="input-pair"),
                       ], className="inputs-block"),

            html.Button("Change Theme", id="toggle-theme-btn", n_clicks=0, className="toggle-btn")

            ], className="sidebar-content")]),




    # ========================================================Main Content ===-------------------------------------------------
    html.H1("Fast Ion Diagnostics Placement Optimizer", className="main-title"),
    html.P("work in progress...", className="wip-label"),

    html.Div([
        html.Div([

            html.Div([
                html.Label("Select Port"),

                dcc.Dropdown(id="nbi-dropdown", options=[{"label": f"NBI_{i}", "value": f"NBI_{i}"} for i in range(1, 9)] + [{"label": f"Gyrotron_{i}", "value": f"Gyrotron_{i}"} for i in range(1, 3)],
           className="custom-dropdown"),

                dcc.Dropdown(id="port-dropdown",options=[],value=None,className="custom-dropdown"),

                html.Label("Diagnostic Type"),

                dcc.Dropdown(id="diagnostic-dropdown", options=[{"label": "FIDA", "value": "FIDA"},{"label": "CTS", "value": "CTS"}],className="custom-dropdown"),],
                    className="custom-dropdown-container"),  

            html.Div([
                dcc.Loading( id="loading-add-port",type="default", children=[
                    html.Button("Add Port and Build", id="add-port-build-btn", n_clicks=0, className="big-action-btn")]),], className="build-buttons"), 
                
                html.Div(className="separator-line"),
                
                html.Div([
                    html.Div([
                        html.Label("Select Config"),
                        dcc.Dropdown(id="config-dropdown",options=[{"label": os.path.basename(path), "value": path} for path in get_config_files()],value=None,
                                     className="custom-dropdown"),], className="tools-block"),
                                     
                        dcc.Upload(id="upload-config", children=html.Button("⬆️ Upload Config", className="big-action-btn"),multiple=False ),
                        html.Div(className="separator-line"),], className="tools-block"),
                        
                    html.Div([
                        # === Pre Calculate ===
                        dcc.Loading(id="loading-precalc",type="default",children=[
                                     html.Button("Pre Calculate", id="build-graph-btn", n_clicks=0, className="big-action-btn"),
                                     html.Div(id="precalc-spinner-output")]) ], className="build-buttons"),


                    html.Div([
                        dcc.Loading(id="build_std",type="default",children=[
                                     html.Button("Build Standard", id='build-standard-btn', n_clicks=0, className="big-action-btn", disabled=True),
                                     html.Div(id="build-standard-output")]),


                        html.Div(className="separator-line") ], className="build-buttons"),
                        
                        
                    html.Div([
                        html.Button("💾 Save Image", id="save-image-btn", n_clicks=0, className="big-action-btn"),
                        dcc.Loading( id="loading-save-results", type="default", children=[
                                html.Button("💾 Save Results", id="save-results-btn", n_clicks=0, className="big-action-btn"),
                                html.Div(id="save-results-output")]),


                        dcc.Upload(id="upload-results",children=html.Button("⬆️ Upload Results", className="big-action-btn"),multiple=False),
                        dcc.Download(id="download-image"),
                        dcc.Download(id="download-results")
                             ], className="tools-block")
        ], className="settings-panel"),

        html.Div([
            html.Img(id='matrix-graph', className='matrix-graph', style={"maxHeight": "700px", 'border-radius': '15px',"overflowY": "auto", "whiteSpace": "pre-wrap", "wordWrap": "break-word"})
            ], className="graph-panel"),

        html.Div([
            html.Pre(id="console-output", children=welcome_message, className="console-panel")
            ], className="console-container"),

    ], id="page-content", className="page-container"),

    dcc.Store(id='theme-store', data={'theme': 'light'}),
    html.Div(id='dummy-output-theme'),
    html.Div(id='dummy-output-sidebar')
])

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
        if (!sidebar) {
            console.log("Sidebar not found!");
            return '';
        }

        if (openClicks == null) openClicks = 0;
        if (closeClicks == null) closeClicks = 0;

        if (openClicks > closeClicks) {
            sidebar.style.width = '250px';
        } else {
            sidebar.style.width = '0';
        }

        return '';
    }
    """,
    Output("dummy-output-sidebar", "children"),
    Input("open-sidebar-btn", "n_clicks"),
    Input("close-sidebar-btn", "n_clicks")
)


# ==========================================================
# 🟢 Python callback: Ports for NBI or Gyrotron
# ==========================================================
@app.callback(
    Output('port-dropdown', 'options'),
    Output('port-dropdown', 'value'),
    Input('nbi-dropdown', 'value'),
    Input('angle-slider', 'value'),
    Input('matrix-size-slider', 'value'),
)
def update_ports_for_nbi(selected_nbi, angle_value, scale_value):
    if selected_nbi is None:
        return [], None

    nbi_index = int(selected_nbi.split("_")[1])
    if selected_nbi.startswith("Gyrotron"):
         nbi_index = nbi_index+8
    nbi_index = nbi_index-1

    _, _, _, valid_port_names = data_instance.port_for_nbi( nbi_index, angle_value, scale_value )

    options = [{'label': name, 'value': name} for name in valid_port_names]
    value = valid_port_names[0] if valid_port_names else None

    return options, value

# ==========================================================
# 🟢 Python callback: Build pre data
# ==========================================================

# here we will prepare all data for start buildins that will have the folling format:
        #====================================================================DATA TYPE ====================================================================================================
        #data_B: [Name NBI; Name Port; Points on NBI; Mag field in this points; Angle between linesight and vec NBI; vec Mag field in points on NBI; angle between vec linesi and magfield]
        #data_B[0]: Name NBI
        #data_B[1]: Name Port
        #data_B[2]: Points on NBI
        #data_B[3]: Mag field in this points
        #data_B[4]: Angle between linesight and vec NBI
        #data_B[5]: vec Mag field in points on NBI
        #data_B[6]: angle between vec linesi and magfield
        #data_B[7]: S
        #data_B[8]: B_max 
        #data_B[9]: J_0
        #data_B[10]: WF
        #==================================================================================================================================================================================


@app.callback(
    Output("calculation-flag", "data", allow_duplicate=True ),
    Input("build-graph-btn", "n_clicks"),
    State("calculation-flag", "data"),
    prevent_initial_call=True
)
def set_calculation_flag(_, is_running):
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True

@app.callback(
    Output("log-list-store", "data", allow_duplicate=True),
    Output("build-graph-btn", "disabled", allow_duplicate=True),
    Output("build-standard-btn", "disabled", allow_duplicate=True ),
    Output("calculation-flag", "data", allow_duplicate=True ), 
    Input("calculation-flag", "data"),
    Input("precalc-flag", "data"),
    State('angle-slider', 'value'),
    State('matrix-size-slider', 'value'),
    State('b0-input', 'value'),
    State('config-dropdown', 'value'), 
    prevent_initial_call=True,
    allow_duplicate=True 
)

def build_graph_with_precalc(is_running, is_running_build, angle, matrix_size, b0, selected_config):
    if not is_running:
        raise dash.exceptions.PreventUpdate
    if is_running_build:
        raise dash.exceptions.PreventUpdate

    os.chdir(BASE_DIR)
    logs = []

    try:
     if selected_config.startswith("CACHED::"):
            filename = selected_config.replace("CACHED::", "")
            all_cached = cache.get(user_cache_key("uploaded_configs")) or {}
            config_data = all_cached.get(filename)
            if config_data is None:
                raise ValueError(f"❗ Cached config '{filename}' not found.")
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode='wb') as tmp_file:
                tmp_file.write(config_data)
                config_path = tmp_file.name
     else:
            config_path = selected_config
     logs.append(f"📁 Config: {os.path.basename(config_path)}")
     logs.append("🚀 Calculating Standart System")
    except Exception as e:
        logs.append(f"📁 Please select a configuration")
        logs.append("❌ Calculation stopped...")

    #main array
    data_full = [[] for _ in range(12)]

    #Standart elements 
    data_full[0] = ['NBI_7','NBI_7', 'NBI_7','NBI_8', 'NBI_8','NBI_8']
    data_full[1] = ['2_1_AEA', '2_1_AEM', '2_1_AET', '2_1_AEA', '2_1_AEM', '2_1_AET']
    data_full[11] = ['FIDA', 'FIDA', 'FIDA', 'FIDA', 'FIDA', 'FIDA']

    try:
        Result_array = data_instance.data_already_input( matrix_size, data_full[1], data_full[0], angle, config_path, b0, data_full[11])
        data_full = Result_array
        logs.append("✅ Done. Points: 3")
    except Exception as e:
        data_full = None
        #logs.append(f"❌ Error: {str(e)}")
    cache.set(user_cache_key("user_result_array"), data_full)
    return logs, False, False, False  




# ==========================================================
# 🟢 Python callback: Build Graph 
# ==========================================================

@app.callback(
    Output("precalc-flag", "data", allow_duplicate=True ),
    Input("build-standard-btn", "n_clicks"),
    State("precalc-flag", "data"),
    prevent_initial_call=True
)
def set_build_flag(n_clicks, is_running):
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True

@app.callback(
    Output("matrix-graph", "src"),
    Output("log-list-store", "data",allow_duplicate=True),
    Output("build-standard-btn", "disabled", allow_duplicate=True ),
    Output("build-standard-output", "children", allow_duplicate=True ),
    Output("precalc-spinner-output", "children", allow_duplicate=True ),
    Output("calculation-flag", "data", allow_duplicate=True),
    Output("precalc-flag", "data", allow_duplicate=True ),
    Input("precalc-flag", "data"),
    State("delta-j0-input", "value"),
    State("delta-s-input", "value"),
    prevent_initial_call=True,
)
def build_matrix_and_update(is_running, delta_J_0, delta_s):
    if not is_running:
        raise dash.exceptions.PreventUpdate

    def replace_none_with_nan(arr_list):
     return [np.array([x if x is not None else np.nan for x in row], dtype=float) for row in arr_list]
    data  = cache.get(user_cache_key("user_result_array"))
    if data is None:
        logs = ["❗ Data not precalculated. Please run Precalculate first."]
        return dash.no_update, logs, False, "",  False,False

    data[9] = replace_none_with_nan(data[9])
    data[10] = replace_none_with_nan(data[10])


    data_wf = data[10][0]
    print(len(data_wf))
    all_results = data
    Name_Ports = data[1]

    logs = []
    logs.append(f"🚀 Starting matrix build for {len(data_wf)} ports")

    num_arrays = len(data_wf)
    fig, axs = plt.subplots(num_arrays, num_arrays, figsize=(10, 10))
    Matr = np.empty((num_arrays, num_arrays), dtype=object)

    args_list = [(i, j, data_wf[i], data_wf[j]) for i in range(num_arrays) for j in range(num_arrays)]
    all_results_list = [all_results] * len(args_list)
    delta_J_0_list = [delta_J_0] * len(args_list)
    delta_s_list = [delta_s] * len(args_list)

    with ProcessPoolExecutor(max_workers=5) as executor:
        results = list(executor.map(calc.compute_matrix, args_list, all_results_list, delta_J_0_list, delta_s_list))

    for i, j, MATRIX, _ in results:
        Matr[i, j] = MATRIX

    min_value, max_value = -1, 3.5

    def plot_subplot(i, j, matrix):
        ax = axs[i, j]
        ax.imshow(matrix, cmap='jet', origin='upper', aspect='auto',
                  vmin=min_value, vmax=max_value)
        ax.set_xticks([])
        ax.set_yticks([])

    with ThreadPoolExecutor(max_workers=5) as executor:
        _ = list(executor.map(
            lambda args: plot_subplot(*args),
            [(i, j, Matr[i, j]) for i in range(num_arrays) for j in range(num_arrays)]
        ))


    for i in range(num_arrays):
        fonts = 6 if num_arrays >= 11 else 9
        selected_nbi = data[0][i]          
        selected_port = Name_Ports[i]  


        if selected_nbi.startswith("NBI"):
            name = "S"
            index_nbi_gyr = selected_nbi.split("_")[1]
        else:
            name = "C"
            index_nbi_gyr = selected_nbi.split("_")[1]


        parts = selected_port.split("_")
        formatted_port = f"{parts[0]}{parts[1]}{parts[2]}"

        label = f"{formatted_port}.{name}{index_nbi_gyr}"

        axs[num_arrays - 1, i].set_xlabel(label, fontsize=fonts)
        axs[i, 0].set_ylabel(label, fontsize=fonts)

    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    plt.colorbar(axs[0, 0].images[0], cax=cax)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)


    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=1500, bbox_inches='tight')
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close(fig)

    img_src = f"data:image/png;base64,{encoded}"
    logs.append("✅ Matrix plot generated successfully.")

    return img_src, logs, False, None, None, False,False

# ==========================================================
# 🟢 Python callback: Update logs
# ==========================================================

@app.callback(
    Output("log-store", "data"),
    Output("console-output", "children"),
    Input("log-list-store", "data"),
    State("log-store", "data"),
    prevent_initial_call=True
)
def update_log_text(new_logs, existing_text):
    if not new_logs:
        raise dash.exceptions.PreventUpdate
    from datetime import datetime
    timestamp = datetime.now().strftime("[%H:%M:%S]")
    lines = "\n".join(f"{timestamp} {line}" for line in new_logs)
    updated = (existing_text or "") + "\n" + lines + "\n"
    return updated, updated



# ==========================================================
# 🟢 Python callback: Save image
# ==========================================================
@app.callback(
    Output("download-image", "data"),
    Input("save-image-btn", "n_clicks"),
    State("matrix-graph", "src"),
    prevent_initial_call=True
)
def download_image(n_clicks, img_src):
    if not img_src:
        raise dash.exceptions.PreventUpdate

    content = img_src.split(',')[1]
    decoded = base64.b64decode(content)

    filename = f"matrix_plot_{datetime.now().strftime('%Y-%m-%d_%H-%M')}.png"
    
    return dcc.send_bytes(decoded, filename=filename)

# ==========================================================
# 🟢 Python callback: Upload new config
# ==========================================================


from dash.exceptions import PreventUpdate

@app.callback(
    Output("config-dropdown", "options", allow_duplicate=True),
    Output("config-dropdown", "value", allow_duplicate=True),
    Output("log-list-store", "data", allow_duplicate=True),
    Input("upload-config", "contents"),
    State("upload-config", "filename"),
    State("log-list-store", "data"),
    prevent_initial_call=True,
)
def handle_config_upload(contents, filename, old_logs):
    if contents is None:
        raise dash.exceptions.PreventUpdate
    
    logs = []

    if not filename.endswith(".txt"):
        logs.append(f"❌ Unsupported file format: {filename}")
        return dash.no_update, dash.no_update, logs

    import base64
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)


    current_configs = cache.get(user_cache_key("uploaded_configs")) or {}

    current_configs[filename] = decoded
    cache.set(user_cache_key("uploaded_configs"), current_configs)

    disk_files = get_config_files()
    options = [{"label": os.path.basename(f), "value": f} for f in disk_files]

    for name in current_configs:
        options.insert(0, {"label": f"(Uploaded) {name}", "value": f"CACHED::{name}"})
    
    logs.append(f"✅ Config uploaded: {filename}")
    return options, f"CACHED::{filename}", logs


# ==========================================================
# 🟢 Python callback: Add new port to grid
# ==========================================================
@app.callback(
    Output("add-port-flag", "data"),
    Input("add-port-build-btn", "n_clicks"),
    State("add-port-flag", "data"),
    prevent_initial_call=True
)
def set_add_port_flag(n_clicks, is_running):
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True

@app.callback(
    Output("log-list-store", "data", allow_duplicate=True),
    Output("add-port-flag", "data", allow_duplicate=True),
    Output("add-port-build-btn", "disabled", allow_duplicate=True),
    Input("add-port-flag", "data"),
    State("nbi-dropdown", "value"),
    State("port-dropdown", "value"),
    State("angle-slider", "value"),
    State("matrix-size-slider", "value"),
    State("config-dropdown", "value"),
    State("b0-input", "value"),
    State("diagnostic-dropdown", "value"),
    State("log-list-store", "data"),
    prevent_initial_call=True
)
def add_port_and_build(is_running, nbi, port, angle, scale, config, B0, diag_type, old_logs):
    if not is_running:
        raise dash.exceptions.PreventUpdate

    logs = []

    data_full = cache.get(user_cache_key("user_result_array"))
    if data_full is None:
        logs.append("❗ Cannot add port: no base data found. Run Pre Calculate first.")
        return logs, False, False



    if not config:
        logs.append("❗ No configuration selected. Please select a config before adding a port.")
        return logs, False, False

    if config.startswith("CACHED::"):
        filename = config.replace("CACHED::", "")
        all_cached = cache.get(user_cache_key("uploaded_configs")) or {}
        config_data = all_cached.get(filename)
        if config_data is None:
            logs.append(f"❌ Cached config '{filename}' not found.")
            return logs, False, False
        import tempfile
        with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode='wb') as tmp_file:
            tmp_file.write(config_data)
            config_path = tmp_file.name
    else:
        config_path = config



    try:
        new_data = data_instance.data_nbi_ports(nbi, port, angle, scale, config_path, B0)
    except Exception as e:
        logs.append(f"❌ Error adding port: {e}")
        return logs, False, False


    for i in range(len(new_data)):
        if i < 11:
            data_full[i].append(new_data[i])

    try:
        i_new = len(data_full[1]) - 1  
        wf_result = calc.process_data_for_point(i_new, data_full, diag_type)
        data_full[10][0].append(wf_result)
    except Exception as e:
        logs.append(f"⚠️ WF not calculated for {port}: {e}")
    print("Updated data_full shape:", [len(x) for x in data_full])


    cache.delete(user_cache_key("user_result_array"))
    cache.set(user_cache_key("user_result_array"), data_full)

    logs.append(f"✅ Added port {port} for {nbi} and calculated WF.")
    return logs, False, False



# ==========================================================
# 🟢 Python callback: Save results
# ==========================================================
import json

@app.callback(
    Output("save-results-flag", "data"),
    Input("save-results-btn", "n_clicks"),
    State("save-results-flag", "data"),
    prevent_initial_call=True
)
def set_save_results_flag(n_clicks, is_running):
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True


@app.callback(
    Output("download-results", "data"),
    Output("log-list-store", "data", allow_duplicate=True),
    Output("save-results-btn", "disabled", allow_duplicate=True),   # ✅ добавили
    Output("save-results-flag", "data", allow_duplicate=True),
    Input("save-results-flag", "data"),
    State("log-list-store", "data"),
    prevent_initial_call=True
)
def save_results(is_running, old_logs):
    if not is_running:
        raise dash.exceptions.PreventUpdate

    logs = []
    data_full = cache.get(user_cache_key("user_result_array"))

    if not data_full or (isinstance(data_full, list) and len(data_full) == 0):
        logs.append("❗ Cannot save results: dataset is empty.")
        return dash.no_update, logs, False, False  # ⬅️ disabled=False, флаг сбрасываем

    def serialize(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, (np.float32, np.float64)):
            return float(obj)
        if isinstance(obj, (np.int32, np.int64)):
            return int(obj)
        return obj

    safe_data = json.dumps(data_full, default=serialize, indent=2)
    filename = f"results_{datetime.now().strftime('%Y-%m-%d_%H-%M')}.json"

    logs.append(f"✅ Results saved as {filename}")
    return dcc.send_string(safe_data, filename=filename), logs, False, False




# ==========================================================
# 🟢 Run the server
# ==========================================================

if __name__ == '__main__':
    HOST = "0.0.0.0"
    PORT = 8050
    local_ip = get_local_ip()
    print("=" * 50)
    print(f" 🚀 Your Dash app is running at:")
    print(f" 👉 Local:   http://127.0.0.1:{PORT}")
    print(f" 👉 Network: http://{local_ip}:{PORT}")
    print("=" * 50)
    app.run(host=HOST, port=PORT, debug=True)
