import dash
from dash import html, dcc, Input, Output, State, clientside_callback
import numpy as np
import matplotlib.pyplot as plt
import io
import base64
import socket
from modules import Data, calculus
import os
import customtkinter as ctk
import threading
import time
from dash import callback_context
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from flask_caching import Cache
from datetime import datetime
import tempfile
from flask import request
import uuid
import base64, tempfile, os
import re  

# ==========================================================
# Server
# ==========================================================

# ---------- Paths ----------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))   # Root directory of the app
CONFIG_DIR = os.path.join(BASE_DIR, "config_files")     # Where configuration files are stored

# ---------- Locks ----------
calc_lock = threading.Lock()  # Prevent concurrent calculation conflicts

# ---------- Network ----------
def get_local_ip():
    """
    Get the local IPv4 address for sharing the app on the local network.
    Falls back to 127.0.0.1 if no external connection is available.
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # This does not actually send data to 8.8.8.8, just used to select a valid network interface
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
    except Exception:
        ip = "127.0.0.1"
    finally:
        s.close()
    return ip

#------------------------- Save Results --------------------
# ===================== HDF5 ragged-safe save/load of full data_full[0..11] =====================
import h5py, numpy as np
from h5py import string_dtype

# Names for readability inside the .h5 file
_H5_NAMES = [
    "nbi_names",              # [0] strings
    "port_names",             # [1] strings
    "nbi_points",             # [2] numerics (ragged ok)
    "mag_field_values",       # [3] numerics
    "angle_linesight_nbi",    # [4] numerics
    "mag_field_vectors",      # [5] numerics
    "angle_linesight_b",      # [6] numerics
    "S_values",               # [7] numerics
    "B_max",                  # [8] numerics
    "J0",                     # [9] numerics (often multi-dim)
    "WF",                     # [10] numerics (Weight Function)
    "diag_types",             # [11] strings
    "config_tags",            # [12] 
    "b0_tags",                # [13] 
]

def _to_numeric_ndarray(x):
    """
    Convert Python structure (list/tuple/ndarray/scalar) to a float32 ndarray,
    replacing None -> np.nan. Works for any nesting/shape.
    """
    arr = np.asarray(x, dtype=object)
    # Replace None with nan, then try to cast to float
    arr = np.vectorize(lambda v: np.nan if v is None else v, otypes=[object])(arr)
    try:
        return arr.astype(np.float32, copy=False)
    except Exception:
        # Last resort: float64
        return arr.astype(np.float64, copy=False)

def _save_string_list(h5, name, lst):
    utf8 = string_dtype(encoding='utf-8')
    data = np.array([str(s) for s in (lst or [])], dtype=object)
    h5.create_dataset(name, data=data, dtype=utf8)

def _save_numeric_list_group(h5, group_name, seq):
    """
    Save a Python list 'seq' of numeric arrays (any shapes) into group/<idx> datasets.
    If seq is None, store empty group.
    """
    grp = h5.create_group(group_name)
    if not seq:
        return
    for i, entry in enumerate(seq):
        arr = _to_numeric_ndarray(entry)
        grp.create_dataset(str(i), data=arr, compression="gzip", shuffle=True, chunks=True)

def _load_string_list(h5, name):
    """Read a string dataset robustly (UTF-8), handling bytes/scalars/arrays."""
    if name not in h5:
        return []
    raw = h5[name][()]

    # normalize to 1D iterable
    arr = np.atleast_1d(raw)
    out = []
    for v in arr:
        if isinstance(v, (bytes, np.bytes_)):
            out.append(v.decode("utf-8", errors="replace"))
        else:
            out.append(str(v))
    return out

def _load_numeric_list_group(h5, group_name):
    """
    Load group of datasets into a Python list (each element is .tolist()).
    Missing group -> [].
    """
    if group_name not in h5:
        return []
    grp = h5[group_name]
    keys = sorted((k for k in grp.keys()), key=lambda x: int(x))  # keep original order
    return [np.array(grp[k][()]).tolist() for k in keys]

def save_data_full_h5(path, data_full):
    """
    Save the FULL 12-slot data_full to HDF5 (ragged-safe).
    Layout:
      - nbi_names, port_names, diag_types: datasets of UTF-8 strings
      - All numeric indices (2..10): groups of datasets (e.g., WF/0, WF/1, ...).
        WF: if your app uses data_full[10] as [ wf_list ], we will flatten to wf_list for storage,
        and restore to [ wf_list ] when loading.
    """
    if not isinstance(data_full, list) or len(data_full) < 14:
        raise ValueError("data_full must be a list of length >= 12")

    with h5py.File(path, "w") as h5:
        h5.attrs["schema_version"] = 3
        # Strings
        _save_string_list(h5, _H5_NAMES[0], data_full[0])   # nbi_names
        _save_string_list(h5, _H5_NAMES[1], data_full[1])   # port_names
        _save_string_list(h5, _H5_NAMES[11], data_full[11]) # diag_types
        _save_string_list(h5, _H5_NAMES[12], data_full[12])
        _save_string_list(h5, _H5_NAMES[13], data_full[13])

        # Numeric groups (2..9 and 10)
        for idx in range(2, 10):
            _save_numeric_list_group(h5, _H5_NAMES[idx], data_full[idx])

        # WF: support both formats: either [ wf_list ] or wf_list itself
        wf_container = data_full[10]
        if isinstance(wf_container, list) and len(wf_container) > 0 and isinstance(wf_container[0], (list, tuple)):
            wf_list = wf_container[0]
            h5.attrs["wf_nested"] = True
        else:
            wf_list = wf_container
            h5.attrs["wf_nested"] = False
        _save_numeric_list_group(h5, _H5_NAMES[10], wf_list)

def load_data_full_h5(path):
    """
    Load FULL 14-slot data_full saved by save_data_full_h5().
    Reconstructs the exact structure, including WF nesting as [ wf_list ].
    For older files (without config_tags/b0_tags), fills with blanks of proper length.
    """
    data_full = [[] for _ in range(14)]  # 0..13
    with h5py.File(path, "r") as h5:

        data_full[0]  = _load_string_list(h5, _H5_NAMES[0])   # nbi_names
        data_full[1]  = _load_string_list(h5, _H5_NAMES[1])   # port_names
        data_full[11] = _load_string_list(h5, _H5_NAMES[11])  # diag_types


        for idx in range(2, 10):
            data_full[idx] = _load_numeric_list_group(h5, _H5_NAMES[idx])


        wf_list = _load_numeric_list_group(h5, _H5_NAMES[10])
        nested = bool(h5.attrs.get("wf_nested", True))
        data_full[10] = [wf_list] if nested else wf_list


        data_full[12] = _load_string_list(h5, _H5_NAMES[12]) if _H5_NAMES[12] in h5 else []
        data_full[13] = _load_string_list(h5, _H5_NAMES[13]) if _H5_NAMES[13] in h5 else []


    ports_n = len(data_full[1]) if isinstance(data_full[1], list) else 0
    data_full[12] = (data_full[12] + [""] * ports_n)[:ports_n]
    data_full[13] = (data_full[13] + [""] * ports_n)[:ports_n]

    return data_full

# ==========================================================
# CONFIGURATION FILES & SESSION KEYS
# ==========================================================

# ---------- List available config files ----------
def get_config_files():
    """
    Return a list of all .txt configuration files in CONFIG_DIR.
    """
    return [
        os.path.join(CONFIG_DIR, f)
        for f in os.listdir(CONFIG_DIR)
        if f.endswith(".txt")
    ]

# ---------- Generate unique cache key per user session ----------
_session_ids = {}

def user_cache_key(base_key="user_result_array"):
    """
    Generate a cache key that is unique for each user session.
    Uses client IP + random UUID (stored for session lifetime).
    """
    try:
        user_ip = request.remote_addr or "anon"
    except RuntimeError:
        # If called outside of request context (e.g., from server), use 'server'
        user_ip = "server"

    if user_ip not in _session_ids:
        _session_ids[user_ip] = str(uuid.uuid4())

    return f"{base_key}_{_session_ids[user_ip]}"

# ==========================================================
# 🟢 PORT LABEL FORMATTING
# ==========================================================

def format_port_labels(port_names, device_names):
    """
    Format port and device names into standardized short labels.

    Example:
        port_names  = ["2_1_AEA", "2_1_AEM"]
        device_names = ["NBI_7", "CTS_1"]
        -> ["21AEA.S7", "21AEM.C1"]

    Rules:
        - Remove underscores from port name (e.g., "2_1_AEA" -> "21AEA")
        - 'S' for NBI (FIDA), 'C' for Gyrotron/CTS
        - Append device index after the letter (e.g., ".S7")
    """
    if not port_names or not device_names:
        return []

    labels = []
    for port, device in zip(port_names, device_names):
        if not port or not device:
            continue

        # Normalize
        port = str(port).strip()
        device = str(device).strip()

        # Port code: remove underscores
        port_code = port.replace("_", "")

        # Device symbol
        dev_symbol = "S" if device.upper().startswith("NBI") else "C"

        # Device index (safe parsing)
        parts = device.split("_", 1)
        dev_index = parts[1] if len(parts) == 2 and parts[1] else parts[0].lstrip("NBI").lstrip("CTS").lstrip("GYROTRON")
        dev_index = dev_index if dev_index else "?"

        labels.append(f"{port_code}.{dev_symbol}{dev_index}")
    return labels


def make_config_tag_from_filename(path_or_name: str) -> str:
    """
    Build a short human-readable config tag from a filename.
    Example: 'EIM_beta=0.txt' -> 'EIM β=0'
    """
    base = os.path.basename(path_or_name or "").strip()
    stem, _ = os.path.splitext(base)
    # First 3 letters as code (fallback to the whole stem if shorter)
    prefix = stem[:3].upper() if stem else "CFG"
    # Extract 'beta=...' (number)
    m = re.search(r'(?i)beta\s*=\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)', stem or "")
    if m:
        return f"{prefix} β={m.group(1)}"
    return prefix

def make_b0_tag(b0) -> str:
    """Format B₀ value for display."""
    try:
        return f"B₀={float(b0):.2f} T"
    except Exception:
        return "B₀=?"


def build_current_config_message(port_names, device_names, diag_types=None, config_tags=None, b0_tags=None):
    """
    Build a formatted 'Current Matrix Configuration' message.
    Adds per-port config tag and B0 tag when available.
    """
    lines = [welcome_message, ""]

    if not port_names or not device_names:
        lines += ["📌 Current Matrix Configuration:", "  ⚠️ No ports loaded. Please Precalculate or Add ports."]
        return "\n".join(lines)

    labels = format_port_labels(port_names, device_names)
    lines.append("📌 Current Matrix Configuration:")

    n = len(labels)
    for i in range(n):
        diag = f" ({str(diag_types[i]).strip()})" if diag_types and i < len(diag_types) and diag_types[i] else ""
        cfg  = f" — {str(config_tags[i]).strip()}" if config_tags and i < len(config_tags) and config_tags[i] else ""
        b0   = f"; {str(b0_tags[i]).strip()}"     if b0_tags and i < len(b0_tags) and b0_tags[i] else ""
        lines.append(f"  • {labels[i]}{diag}{cfg}{b0}")

    return "\n".join(lines)

# ==========================================================
# CONSTANTS
# ==========================================================

welcome_message = (
    "Explanation of the naming convention:\n"
    "For '21AEM.S8':\n"
    "  - 21:  Module 2, Submodule 1\n"
    "  - AEM: Type of the port\n"
    "  - S:   NBI (FIDA)\n"
    "  - C:   Gyrotron (CTS)\n"
    "  - 8:   NBI or Gyrotron number\n"
    "-------------------------------------"
)

# ==========================================================
# INITIALIZE 
# ==========================================================

# Create Dash app instance
app = dash.Dash(
    __name__,
    suppress_callback_exceptions=True,  # Allow callbacks for components not in initial layout
    prevent_initial_callbacks=True      # Prevent firing callbacks on initial page load
)
app.enable_duplicate_callbacks = True  # Allow multiple callbacks to have same Output

# Configure Flask-Caching (using filesystem backend)
cache = Cache(
    app.server,
    config={
        'CACHE_TYPE': 'filesystem',
        'CACHE_DIR': 'cache-dir',
        'CACHE_DEFAULT_TIMEOUT': 3600  # seconds
    }
)
# Initialize data/calculation modules
data_instance = Data()
calc = calculus()

# ==========================================================
# 🟢 APP LAYOUT
# ==========================================================
app.layout = html.Div([

    # ====================== STATE STORAGE ======================
    dcc.Store(id="calculation-flag", data=False),
    dcc.Store(id="log-list-store"),
    dcc.Store(id="log-store", data=""),
    dcc.Store(id="precalc-flag", data=False),
    dcc.Store(id="standard-flag", data=False),
    dcc.Store(id="add-port-flag", data=False),
    dcc.Store(id="save-results-flag", data=False),
    dcc.Store(id="cache-update-trigger", data=0),

    # Interval for initial header rendering
    dcc.Interval(id="header-init-interval", interval=1, n_intervals=0, max_intervals=1),

    # ====================== SIDEBAR TOGGLE BUTTON ======================
    html.Span("☰", id="open-sidebar-btn", className="open-btn"),

    # ====================== SIDEBAR ======================
    html.Div(
        id="sidebar",
        className="sidebar",
        children=[
            html.Div(
                className="sidebar-content",
                children=[

                    # Header
                    html.Div(
                        className="sidebar-header",
                        children=[
                            html.H2([html.Span("⚙️"), " Settings"], className="sidebar-title-inline"),
                            html.Button("☰", id="close-sidebar-btn", className="close-btn"),
                        ],
                    ),

                    # Sliders
                    html.Div(
                        className="slider-block",
                        children=[
                            html.Label("Angle", className="label"),
                            dcc.Slider(
                                id="angle-slider",
                                min=0, max=90, step=1, value=90,
                                tooltip={"always_visible": False, "placement": "bottom"},
                                marks={0: "0°", 90: "90°"},
                            ),
                            html.Label("Matrix Size", className="label"),
                            dcc.Slider(
                                id="matrix-size-slider",
                                min=10, max=100, step=5, value=10,
                                tooltip={"always_visible": False, "placement": "bottom"},
                                marks={10: "10", 50: "50", 100: "100"},
                            ),
                        ],
                    ),

                    # Inputs
                    html.Div(
                        className="inputs-block",
                        children=[
                            html.Div(
                                className="input-pair",
                                children=[
                                    html.Label("Δs =", className="inline-label"),
                                    dcc.Input(id="delta-s-input", type="number", value=0.05, step=0.01, className="sidebar-input"),
                                ],
                            ),
                            html.Div(
                                className="input-pair",
                                children=[
                                    html.Label("ΔJ₀ =", className="inline-label"),
                                    dcc.Input(id="delta-j0-input", type="number", value=0.005, step=0.001, className="sidebar-input"),
                                ],
                            ),
                            html.Div(
                                className="input-pair",
                                children=[
                                    html.Label("B₀ =", className="inline-label"),
                                    dcc.Input(id="b0-input", type="number", value=2.52, step=0.01, className="sidebar-input"),
                                ],
                            ),
                        ],
                    ),

                    # Theme switch
                    html.Button("Change Theme", id="toggle-theme-btn", n_clicks=0, className="toggle-btn"),
                ],
            ),
        ],
    ),

    # ====================== MAIN CONTENT ======================
    html.H1("Fast Ion Diagnostics Placement Optimizer", className="main-title"),
    html.P("work in progress...", className="wip-label"),

    html.Div(
        id="page-content",
        className="page-container",
        children=[

            # ---------- LEFT SETTINGS PANEL ----------
            html.Div(
                className="settings-panel",
                children=[

                    # Port + Diagnostic
                    html.Div(
                        className="custom-dropdown-container",
                        children=[
                            html.Label("Select Port"),
                            dcc.Dropdown(
                                id="nbi-dropdown",
                                options=(
                                    [{"label": f"NBI_{i}", "value": f"NBI_{i}"} for i in range(1, 9)]
                                    + [{"label": f"Gyrotron_{i}", "value": f"Gyrotron_{i}"} for i in range(1, 3)]
                                ),
                                className="custom-dropdown",
                            ),
                            dcc.Dropdown(id="port-dropdown", options=[], value=None, className="custom-dropdown"),
                            html.Label("Diagnostic Type"),
                            dcc.Dropdown(
                                id="diagnostic-dropdown",
                                options=[{"label": "FIDA", "value": "FIDA"}, {"label": "CTS", "value": "CTS"}],
                                className="custom-dropdown",
                            ),
                        ],
                    ),

                    # Add Port
                    html.Div(
                        className="build-buttons",
                        children=[
                            dcc.Loading(
                                id="loading-add-port",
                                type="default",
                                children=[html.Button("Add Port", id="add-port-build-btn", n_clicks=0, className="big-action-btn")],
                            ),
                        ],
                    ),

                    html.Div(className="separator-line"),

                    # Config selection
                    html.Div(
                        className="tools-block",
                        children=[
                            html.Label("Select Config"),
                            dcc.Dropdown(
                                id="config-dropdown",
                                options=[{"label": os.path.basename(path), "value": path} for path in get_config_files()],
                                value=None,
                                className="custom-dropdown",
                            ),
                        ],
                    ),
                    dcc.Upload(id="upload-config", children=html.Button("⬆️ Upload Config", className="big-action-btn"), multiple=False),

                    html.Div(className="separator-line"),

                    # --- Matrix pipeline card ---
                    
                    html.Div(
                        className="action-card",
                        children=[
                            html.Div(className="action-card__title", children=["Matrix Pipeline"]),


                            html.Div(
                                className="action-row",
                                children=[
                                    dcc.Loading(
                                        id="loading-precalc",
                                        type="default",
                                        children=[
                                            html.Button(
                                                "Precalculate",
                                                id="build-graph-btn",
                                                n_clicks=0,
                                                className="big-action-btn",
                                                title="Prepare base arrays and store them in cache"
                                            ),
                                            html.Div(id="precalc-spinner-output")
                                        ],
                                    ),

                                    dcc.Loading(
                                        id="build_std",
                                        type="default",
                                        children=[
                                            html.Button(
                                                "Build Grid",
                                                id="build-standard-btn",
                                                n_clicks=0,
                                                className="big-action-btn",
                                                disabled=True,
                                                title="Generate the graph matrix after Precalculate"
                                            ),
                                            html.Div(id="build-standard-output")
                                        ],
                                    ),
                                ],
                            ),
                        ],
                    ),



                    html.Div(className="separator-line"),
                    html.Div(
                        className="tools-block",
                        children=[
                            html.Button("💾 Save Image", id="save-image-btn", n_clicks=0, className="big-action-btn"),

                            # Save Results already has Loading 
                            dcc.Loading(
                                id="loading-save-results",
                                type="default",
                                children=[
                                    html.Button("💾 Save Results", id="save-results-btn", n_clicks=0, className="big-action-btn"),
                                    html.Div(id="save-results-output"),       # <- placeholder for spinner
                                ],
                            ),

                            # ⬇️ NEW: wrap Upload Results with Loading just like Precalc
                            dcc.Loading(
                                id="loading-upload-results",
                                type="default",
                                children=[
                                    dcc.Upload(
                                        id="upload-results",
                                        children=html.Button("⬆️ Upload Results", className="big-action-btn"),
                                        multiple=False,
                                        accept=".h5,.hdf5"  
                                    ),
                                    html.Div(id="upload-results-output"),      # <- placeholder for spinner
                                ],
                            ),

                            dcc.Download(id="download-image"),
                            dcc.Download(id="download-results"),
                        ],
                    ),
                ],
            ),

            # ---------- GRAPH PANEL ----------
            html.Div(
                className="graph-panel",
                children=[
                    html.Img(
                        id="matrix-graph",
                        className="matrix-graph",
                        style={
                            "maxHeight": "700px",
                            "borderRadius": "15px",
                            "overflowY": "auto",
                            "whiteSpace": "pre-wrap",
                            "wordWrap": "break-word",
                        },
                    )
                ],
            ),

            # ---------- CONSOLE ----------
            html.Div(
                className="console-container",
                children=[
                    html.Div(id="console-header", className="console-block console-header"),
                    html.Div(id="console-output", className="console-block console-log"),
                ],
            ),
        ],
    ),

    # Theme store & dummies for clientside callbacks
    dcc.Store(id="theme-store", data={"theme": "light"}),
    html.Div(id="dummy-output-theme"),
    html.Div(id="dummy-output-sidebar"),
])

# ==========================================================
# Clientside callback: Apply theme (light or dark)
# ==========================================================
clientside_callback(
    """
        function(data, openClicks, closeClicks) {
        // --- Theme switcher 
        const theme = (data && data.theme) ? data.theme : 'light';
        document.body.classList.remove('light-theme', 'dark-theme');
        document.body.classList.add(`${theme}-theme`);

        const sidebar = document.getElementById('sidebar');
        const pageContent = document.getElementById('page-content');

        if (sidebar) {
            sidebar.classList.remove('light-sidebar', 'dark-sidebar');
            sidebar.classList.add(theme === 'light' ? 'light-sidebar' : 'dark-sidebar');
        }

        // --- Sidebar shift logic
        openClicks = openClicks || 0;
        closeClicks = closeClicks || 0;

        const isOpen = openClicks > closeClicks;

        if (sidebar) sidebar.style.width = isOpen ? '250px' : '0';
        if (pageContent) pageContent.style.marginLeft = isOpen ? '250px' : '0';

        return ['', ''];
        }
    """,
    Output('dummy-output-theme', 'children'),
    Input('theme-store', 'data'),
    Input('open-sidebar-btn', 'n_clicks'),
    Input('close-sidebar-btn', 'n_clicks')
)


# ==========================================================
# Python callback: Toggle theme and update store
# ==========================================================
@app.callback(
    Output('theme-store', 'data'),
    Input('toggle-theme-btn', 'n_clicks'),
    State('theme-store', 'data'),
    prevent_initial_call=True
)
def toggle_theme(_, data):
    # Get current theme or set to 'light' if not defined
    current_theme = (data or {}).get('theme', 'light')
    
    # Switch between light and dark theme
    new_theme = 'dark' if current_theme == 'light' else 'light'
    
    # Return updated theme data to be stored
    return {'theme': new_theme}

# ==========================================================
# Clientside callback: Open / Close Sidebar
# ==========================================================

clientside_callback(
    """
    function(openClicks, closeClicks) {
        // Get the sidebar element
        const sidebar = document.getElementById('sidebar');
        if (!sidebar) {
            console.warn("Sidebar element not found!");
            return '';
        }

        // Default click counts to 0 if undefined
        if (openClicks == null) openClicks = 0;
        if (closeClicks == null) closeClicks = 0;

        // Toggle sidebar visibility based on clicks
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
# Python callback: Update available ports for selected NBI or Gyrotron
# ==========================================================
@app.callback(
    Output('port-dropdown', 'options'),
    Output('port-dropdown', 'value'),
    Input('nbi-dropdown', 'value'),
    Input('angle-slider', 'value'),
    Input('matrix-size-slider', 'value'),
)
def update_ports_for_nbi(selected_nbi, angle_value, scale_value):
    """
    Updates the 'Port' dropdown based on the selected NBI or Gyrotron.

    Args:
        selected_nbi (str): Selected device, e.g., 'NBI_7' or 'Gyrotron_1'.
        angle_value (int): Current angle slider value.
        scale_value (int): Current matrix size slider value.

    Returns:
        tuple:
            - list[dict]: Dropdown options for available ports.
            - str|None: Default selected port.
    """
    # No selection → clear dropdown
    if selected_nbi is None:
        return [], None

    # Determine index for NBI or Gyrotron
    nbi_index = int(selected_nbi.split("_")[1])
    if selected_nbi.startswith("Gyrotron"):
        nbi_index += 8  # Gyrotrons follow after NBI indices
    nbi_index -= 1  # Zero-based index for internal use

    # Retrieve list of valid ports
    _, _, _, valid_port_names = data_instance.port_for_nbi(nbi_index, angle_value, scale_value)

    # Prepare dropdown data
    options = [{'label': name, 'value': name} for name in valid_port_names]
    value = valid_port_names[0] if valid_port_names else None

    return options, value


# ==========================================================
# Python callback: Build pre data
# ==========================================================

# ==========================================================
# Data Structure Reference: data_B / data_full
# ==========================================================
# The pre-calculation step prepares a 2D list (array) containing
# base data for the diagnostic system configuration.
#
# Index mapping:
# ----------------------------------------------------------
# [0]  nbi_names             → List of NBI or CTS
# [1]  port_names            → List of port names 
# [2]  nbi_points            → Points along NBI lines of sight
# [3]  mag_field_values      → Magnetic field magnitude at those points
# [4]  angle_linesight_nbi   → Angle between line-of-sight vector and NBI vector
# [5]  mag_field_vectors     → Magnetic field vector at those points
# [6]  angle_linesight_b     → Angle between line-of-sight vector and magnetic field vector
# [7]  S_values              → 'S' parameter values
# [8]  B_max                 → Maximum magnetic field 
# [9]  J_0                   → J₀ values for each point
# [10] WF                    → Weight Function
# [11] diag_types            → Diagnostic type 
# [12] config_tags          → Per-port config tag (e.g., "EIM β=0")
# [13] b0_tags              → Per-port B₀ tag (e.g., "B₀=2.52 T")
# ----------------------------------------------------------

@app.callback(
    Output("calculation-flag", "data", allow_duplicate=True),
    Input("build-graph-btn", "n_clicks"),
    State("calculation-flag", "data"),
    prevent_initial_call=True
)
def set_calculation_flag(_, is_running):
    """
    Trigger the calculation flag when 'Precalculate' button is pressed,
    but only if a calculation is not already running.
    """
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True


@app.callback(
    Output("log-list-store", "data", allow_duplicate=True),
    Output("build-graph-btn", "disabled", allow_duplicate=True),
    Output("build-standard-btn", "disabled", allow_duplicate=True),
    Output("calculation-flag", "data", allow_duplicate=True),
    Output("console-header", "children", allow_duplicate=True),
    Input("calculation-flag", "data"),
    Input("precalc-flag", "data"),
    State("angle-slider", "value"),
    State("matrix-size-slider", "value"),
    State("b0-input", "value"),
    State("config-dropdown", "value"),
    prevent_initial_call=True
)
def build_graph_with_precalc(is_running, is_running_build, angle, matrix_size, b0, selected_config):
    """
    Run the precalculation step:
    - Load selected configuration (from disk or cache)
    - Build the initial data arrays
    - Store them in cache for later 'Build Grid'
    """
    if not is_running or is_running_build:
        raise dash.exceptions.PreventUpdate

    logs = []
    config_path = None

    # --- Load configuration ---
    try:
        if not selected_config:
            raise ValueError("No configuration selected")

        if selected_config.startswith("CACHED::"):
            filename = selected_config.replace("CACHED::", "")
            all_cached = cache.get(user_cache_key("uploaded_configs")) or {}
            config_data = all_cached.get(filename)
            if config_data is None:
                raise ValueError(f"Cached config '{filename}' not found")
            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode="wb") as tmp_file:
                tmp_file.write(config_data)
                config_path = tmp_file.name
        else:
            config_path = selected_config

        logs.append(f"📁 Config: {os.path.basename(config_path)}")
    except Exception as e:
        logs.append("📁 Please select a configuration")
        logs.append(f"❌ Calculation stopped: {e}")
        lines = build_current_config_message(None, None, diag_types=None)
        return logs, False, False, False, lines

    # --- Initialize base data array ---
    data_full = [[] for _ in range(14)]
    data_full[0] = ["NBI_7", "NBI_7", "NBI_7", "NBI_8", "NBI_8", "NBI_8"]
    data_full[1] = ["2_1_AEA", "2_1_AEM", "2_1_AET", "2_1_AEA", "2_1_AEM", "2_1_AET"]
    data_full[11] = ["FIDA"] * 6

    # --- Run main data preparation ---
    try:
        result_array = data_instance.data_already_input(
            matrix_size, data_full[1], data_full[0], angle, config_path, b0, data_full[11]
        )
        logs.append("✅ Precalculate complete. Points: 3")
        result_array.append(data_full[11])
        # choose a representative name for config tag
        cfg_name_for_tag = filename if (selected_config.startswith("CACHED::")) else os.path.basename(config_path)
        cfg_tag = make_config_tag_from_filename(cfg_name_for_tag)
        b0_tag  = make_b0_tag(b0)

        num_ports = len(result_array[1]) if (isinstance(result_array, list) and len(result_array) > 1) else 0
        config_tags = [cfg_tag] * num_ports
        b0_tags     = [b0_tag]  * num_ports

        # indices [12] and [13]
        result_array.append(config_tags)
        result_array.append(b0_tags)

        # header with config & B0
        lines = build_current_config_message(
            result_array[1], result_array[0],
            diag_types=result_array[11],
            config_tags=result_array[12],
            b0_tags=result_array[13],
        )
    except Exception as e:
        result_array = None
        logs.append(f"❌ Error during precalculate: {e}")
        lines = build_current_config_message(None, None, diag_types=None)

    # --- Store results in cache ---
    cache.set(user_cache_key("user_result_array"), result_array)

    return logs, False, False, False, lines


# ==========================================================
# Python callback: Build Graph 
# ==========================================================

@app.callback(
    Output("precalc-flag", "data", allow_duplicate=True),
    Input("build-standard-btn", "n_clicks"),
    State("precalc-flag", "data"),
    prevent_initial_call=True
)
def set_build_flag(_, is_running):
    """
    Arm the 'Build Grid' step when the button is pressed,
    but only if there isn't an active build already.
    """
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True

@app.callback(
    Output("matrix-graph", "src"),
    Output("log-list-store", "data", allow_duplicate=True),
    Output("build-standard-btn", "disabled", allow_duplicate=True),
    Output("build-standard-output", "children", allow_duplicate=True),
    Output("precalc-spinner-output", "children", allow_duplicate=True),
    Output("calculation-flag", "data", allow_duplicate=True),
    Output("precalc-flag", "data", allow_duplicate=True),
    Output("cache-update-trigger", "data", allow_duplicate=True),
    Input("precalc-flag", "data"),
    State("delta-j0-input", "value"),
    State("delta-s-input", "value"),
    State("calculation-flag", "data"),
    prevent_initial_call=True,
)
def build_matrix_and_update(is_running, delta_j0, delta_s, trgr):
    """
    Build the WF cross-matrix plot:
      1) Read precalculated data from cache
      2) Clean None → NaN for numeric arrays
      3) Compute the full pairwise matrix in parallel
      4) Render a composite image and return it as base64
      5) Log the current configuration using helper builders
    """
    if not is_running:
        raise dash.exceptions.PreventUpdate

    # --- Load data from cache ---
    data = cache.get(user_cache_key("user_result_array"))
    if data is None:
        logs = ["❗ Data not precalculated. Please run Precalculate first."]
        return dash.no_update, logs, False, None, None, False, False, trgr

    # --- Helpers: replace None with NaN for numeric safety ---
    def replace_none_with_nan(arr_list):
        return [
            np.array([x if x is not None else np.nan for x in row], dtype=float)
            for row in arr_list
        ]

    # Normalize J0 and WF containers
    data[9] = replace_none_with_nan(data[9])
    if data[10] and isinstance(data[10], list) and len(data[10]) > 0:
        wf_list = data[10][0]
        data[10][0] = [
            (np.where(np.equal(wf, None), np.nan, wf).astype(float, copy=False)
            if isinstance(wf, np.ndarray)
            else np.array([np.nan if v is None else v for v in wf], dtype=float))
            for wf in wf_list
        ]

    # --- Extract inputs for matrix build ---
    wf_sets = data[10][0]          # list of WF arrays per port
    port_names = data[1]           # e.g., ["2_1_AEA", ...]
    device_names = data[0]         # e.g., ["NBI_7", ...]
    num_arrays = len(wf_sets)

    logs = []

    # --- Prepare figure and containers ---
    fig, axs = plt.subplots(num_arrays, num_arrays, figsize=(10, 10))
    matrix_cells = np.empty((num_arrays, num_arrays), dtype=object)

    # Build argument lists for the parallel compute
    args_list = [(i, j, wf_sets[i], wf_sets[j]) for i in range(num_arrays) for j in range(num_arrays)]
    all_results_list = [data] * len(args_list)
    delta_j0_list = [delta_j0] * len(args_list)
    delta_s_list = [delta_s] * len(args_list)

    # --- Compute matrix in parallel (process pool) ---
    with ProcessPoolExecutor(max_workers=5) as executor:
        results = list(executor.map(calc.compute_matrix, args_list, all_results_list, delta_j0_list, delta_s_list))

    for i, j, matrix, _ in results:
        matrix_cells[i, j] = matrix

    # --- Plot heatmap grid (thread pool for snappy rendering) ---
    vmin, vmax = -1, 3.5

    def plot_subplot(i, j, mat):
        ax = axs[i, j]
        im = ax.imshow(mat, cmap="jet", origin="upper", aspect="auto", vmin=vmin, vmax=vmax)
        ax.set_xticks([])
        ax.set_yticks([])
        return im

    with ThreadPoolExecutor(max_workers=5) as executor:
        ims = list(executor.map(
            lambda args: plot_subplot(*args),
            [(i, j, matrix_cells[i, j]) for i in range(num_arrays) for j in range(num_arrays)]
        ))

    # --- Axis labels using our shared formatter ---
    # Build compact labels like "21AEM.S8" consistently across app
    labels = format_port_labels(port_names, device_names)
    tick_font = 6 if num_arrays >= 11 else 9

    for i, lbl in enumerate(labels):
        axs[num_arrays - 1, i].set_xlabel(lbl, fontsize=tick_font)
        axs[i, 0].set_ylabel(lbl, fontsize=tick_font)

    # --- Colorbar & layout ---
    cax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    plt.colorbar(axs[0, 0].images[0], cax=cax)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)

    # --- Encode figure to base64 data URI ---
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=1500, bbox_inches="tight")
    buf.seek(0)
    encoded = base64.b64encode(buf.read()).decode("utf-8")
    buf.close()
    plt.close(fig)

    img_src = f"data:image/png;base64,{encoded}"

    # Logs
    logs.append("✅ Matrix plot generated successfully.")
    logs.append("📌 Current Matrix Configuration:")
    logs.extend(f"  • {lbl}" for lbl in labels)

    # Reset both flags; bump cache trigger to refresh any dependents
    return img_src, logs, False, None, None, False, False, trgr + 1


# ==========================================================
# Python callback: Initialize console header
# ==========================================================
@app.callback(
    Output("console-header", "children", allow_duplicate=True),
    Input("header-init-interval", "n_intervals"),
)
def init_header(_):
    """
    Populate the console header with the static welcome message
    when the app starts. This runs only once via header-init-interval.
    """
    return [welcome_message, ""]

# ==========================================================
# 🟢 Python callback: Update logs
# ==========================================================

@app.callback(
    Output("console-output", "children"),
    Output("log-store", "data"),
    Input("log-list-store", "data"),
    State("log-store", "data"),
)
def update_log_output(new_logs, existing_text):
    """
    Prepend new log lines (with timestamps) to the existing console text.
    - Accepts either list[str] or str for `new_logs`
    - Removes the static welcome message from the rolling log once
    - Keeps the log size under control by trimming to the last N lines
    """
    from datetime import datetime

    if not new_logs:
        raise dash.exceptions.PreventUpdate

    # --- Normalize inputs ---
    if isinstance(new_logs, str):
        # Split multiline string into lines
        new_logs = [line for line in new_logs.splitlines() if line.strip()]
    else:
        # Ensure list of non-empty strings
        new_logs = [str(line).strip() for line in new_logs if str(line).strip()]

    if not new_logs:
        raise dash.exceptions.PreventUpdate

    existing_text = existing_text or ""

    # --- Remove welcome block if it slipped into the rolling log previously ---
    if welcome_message.strip() in existing_text:
        existing_text = existing_text.replace(welcome_message.strip(), "").strip()

    # --- Timestamped prefix for each new line ---
    ts = datetime.now().strftime("[%H:%M:%S]")
    new_chunk = "\n".join(f"{ts} {line}" for line in new_logs)

    # --- Prepend new lines to old content ---
    updated = (new_chunk + ("\n\n" + existing_text if existing_text else "")).strip()

    # --- (Optional) Trim the console to last N lines to avoid runaway growth ---
    MAX_LINES = 200  # tune as you like
    lines = updated.splitlines()
    if len(lines) > MAX_LINES:
        updated = "\n".join(lines[:MAX_LINES])  # we keep the newest first (prepended)

    return updated, updated

# ==========================================================
# Python callback: Save image
# ==========================================================
@app.callback(
    Output("download-image", "data"),
    Input("save-image-btn", "n_clicks"),
    State("matrix-graph", "src"),
    prevent_initial_call=True
)
def download_image(_, img_src):
    """
    Download the rendered matrix image displayed in <img id="matrix-graph">.
    Expects a data URI like: data:image/png;base64,AAAA...
    """
    if not img_src:
        raise dash.exceptions.PreventUpdate

    # Ensure this is a data URL with an image payload
    if not img_src.startswith("data:image/"):
        # Nothing to download if src is not an inline image
        raise dash.exceptions.PreventUpdate

    # Parse mimetype and base64 payload
    header, b64data = img_src.split(",", 1)
    # header example: "data:image/png;base64"
    mime = header[5:].split(";")[0]  # -> "image/png"

    # Derive file extension from MIME (fallback to .png)
    ext = {
        "image/png": "png",
        "image/jpeg": "jpg",
        "image/webp": "webp",
        "image/svg+xml": "svg"
    }.get(mime, "png")

    # Decode base64 safely
    try:
        decoded = base64.b64decode(b64data, validate=True)
    except Exception:
        # If payload is invalid, do not trigger a broken download
        raise dash.exceptions.PreventUpdate

    # Timestamped filename
    filename = f"matrix_plot_{datetime.now().strftime('%Y-%m-%d_%H-%M')}.{ext}"

    # Return bytes to trigger browser download
    return dcc.send_bytes(decoded, filename=filename)


# ==========================================================
# Python callback: Upload new config
# ==========================================================

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
    """
    Handle a single .txt config upload:
      - Validate extension and data URL
      - Decode base64 payload
      - Store bytes in per-session cache under a unique key
      - Refresh the dropdown options (uploaded first, then on-disk files)
      - Select the just-uploaded file
    """
    if not contents or not filename:
        raise dash.exceptions.PreventUpdate

    logs = []

    # --- Validate extension (case-insensitive) ---
    if not str(filename).lower().endswith(".txt"):
        logs.append(f"❌ Unsupported file format: {filename}")
        return dash.no_update, dash.no_update, logs

    # --- Parse and decode data URL safely ---
    try:
        if "," not in contents:
            raise ValueError("Malformed data URL")
        header, content_b64 = contents.split(",", 1)
        # Optional: ensure it's text (not strictly required if you trust the UI)
        if not header.startswith("data:"):
            raise ValueError("Invalid content header")
        decoded = base64.b64decode(content_b64, validate=True)
    except Exception as e:
        logs.append(f"❌ Failed to decode uploaded file: {e}")
        return dash.no_update, dash.no_update, logs

    # --- Load current per-session uploaded configs from cache ---
    current_configs = cache.get(user_cache_key("uploaded_configs")) or {}

    # --- Sanitize and uniquify filename within this session cache ---
    safe_name = os.path.basename(filename).strip() or "config.txt"
    base, ext = os.path.splitext(safe_name)
    candidate = safe_name
    suffix = 1
    while candidate in current_configs:
        candidate = f"{base} ({suffix}){ext}"
        suffix += 1
    safe_name = candidate

    # --- Save bytes in cache ---
    current_configs[safe_name] = decoded
    cache.set(user_cache_key("uploaded_configs"), current_configs)

    # --- Build dropdown options: uploaded first (sorted), then disk files (sorted) ---
    disk_files = sorted(get_config_files(), key=lambda p: os.path.basename(p).lower())
    options = []

    # Uploaded (session) configs
    for name in sorted(current_configs.keys(), key=str.lower):
        options.append({"label": f"(Uploaded) {name}", "value": f"CACHED::{name}"})

    # On-disk configs
    options.extend({"label": os.path.basename(p), "value": p} for p in disk_files)

    logs.append(f"✅ Config uploaded: {safe_name}")
    return options, f"CACHED::{safe_name}", logs


# ==========================================================
# Python callback: Add new port to grid
# ==========================================================
@app.callback(
    Output("add-port-flag", "data"),
    Input("add-port-build-btn", "n_clicks"),
    State("add-port-flag", "data"),
    prevent_initial_call=True
)
def set_add_port_flag(n_clicks, is_running):
    """Arm the 'Add Port' action unless it is already running."""
    if is_running:
        raise dash.exceptions.PreventUpdate
    return True


@app.callback(
    Output("log-list-store", "data", allow_duplicate=True),
    Output("add-port-flag", "data", allow_duplicate=True),
    Output("add-port-build-btn", "disabled", allow_duplicate=True),
    Output("build-standard-output", "children", allow_duplicate=True),
    Output("precalc-spinner-output", "children", allow_duplicate=True),
    Output("console-header", "children", allow_duplicate=True),  
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
def add_port_and_build(is_running, nbi, port, angle, scale, config, b0, diag_type, old_logs):
    """
    Add a new (NBI/Gyrotron, Port) to the grid:
      1) Validate inputs and resolve config path (disk or session cache)
      2) Compute new data for the selected port and append to `data_full`
      3) Compute WF for the new point and append to data_full[10][0]
      4) Save updated `data_full` back to cache
      5) Update logs and refresh console header (Current Matrix Configuration)
    """
    if not is_running:
        raise dash.exceptions.PreventUpdate

    logs = []

    # --- Load base data prepared by Precalculate ---
    data_full = cache.get(user_cache_key("user_result_array"))
    if data_full is None:
        logs.append("❗ Cannot add port: no base data found. Run Precalculate first.")
        # keep UI enabled; do not change header
        return logs, False, False, None, None, dash.no_update

    # --- Basic input validation ---
    if not config:
        logs.append("❗ No configuration selected. Please select a config before adding a port.")
        return logs, False, False, None, None, dash.no_update
    if not nbi or not port:
        logs.append("❗ Please select both device and port.")
        return logs, False, False, None, None, dash.no_update

    # --- Resolve config path (disk or session-cache) ---
    try:
        if str(config).startswith("CACHED::"):
            filename = config.replace("CACHED::", "", 1)
            all_cached = cache.get(user_cache_key("uploaded_configs")) or {}
            config_data = all_cached.get(filename)
            if config_data is None:
                logs.append(f"❌ Cached config '{filename}' not found.")
                return logs, False, False, None, None, dash.no_update

            with tempfile.NamedTemporaryFile(delete=False, suffix=".txt", mode="wb") as tmp_file:
                tmp_file.write(config_data)
                config_path = tmp_file.name
        else:
            config_path = config
    except Exception as e:
        logs.append(f"❌ Failed to prepare configuration: {e}")
        return logs, False, False, None, None, dash.no_update

    # --- Compute data for the new port and append to data_full ---
    try:
        new_data = data_instance.data_nbi_ports(nbi, port, angle, scale, config_path, b0)
    except Exception as e:
        logs.append(f"❌ Error adding port: {e}")
        return logs, False, False, None, None, dash.no_update

    # Append elements 0..10 (index 11 is diag types list)
    for i in range(min(len(new_data), 11)):
        data_full[i].append(new_data[i])

    # Ensure diag type list exists & append diagnostic type for the new entry
    if len(data_full) < 14:
        data_full.extend([[] for _ in range(14 - len(data_full))])
    if not data_full[11]:
        data_full[11] = []
    data_full[11].append(diag_type or "FIDA")

    # --- Compute WF for the new point and append to data_full[10][0] ---
    try:
        i_new = len(data_full[1]) - 1  # index of the newly added port
        wf_result = calc.process_data_for_point(i_new, data_full, diag_type)
        # data_full[10] is expected to be a list where [0] is the WF list
        if not data_full[10]:
            data_full[10] = [[]]
        if not isinstance(data_full[10][0], list):
            data_full[10][0] = list(data_full[10][0])
        data_full[10][0].append(wf_result)
        # Append per-port config & B0 tags
        cfg_name_for_tag = config.replace("CACHED::", "", 1) if str(config).startswith("CACHED::") else os.path.basename(config)
        data_full[12].append(make_config_tag_from_filename(cfg_name_for_tag))
        data_full[13].append(make_b0_tag(b0))
    except Exception as e:
        logs.append(f"⚠️ WF not calculated for {port}: {e}")

    # --- Persist updated data to cache ---
    cache.set(user_cache_key("user_result_array"), data_full)

    # --- Logs & header update ---
    logs.append(f"✅ Added port {port} for {nbi} and calculated WF.")

    # Refresh header using the shared builder (includes welcome + nicely formatted labels)
    try:
        header_lines = build_current_config_message(
            port_names=data_full[1],
            device_names=data_full[0],
            diag_types=data_full[11],
            config_tags=data_full[12],
            b0_tags=data_full[13],
        )
    except Exception:
        # Fallback: keep previous header unchanged if formatter fails
        header_lines = dash.no_update

    # Disable the flag and button, clear spinners; update header
    return logs, False, False, None, None, header_lines




# ==========================================================
# 🟢 Python callback: Save results to HDF5
# ==========================================================

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
    Output("download-results", "data", allow_duplicate=True),
    Output("log-list-store", "data", allow_duplicate=True),
    Output("save-results-btn", "disabled", allow_duplicate=True),
    Output("save-results-flag", "data", allow_duplicate=True),
    Input("save-results-flag", "data"),
    State("log-list-store", "data"),
    prevent_initial_call=True
)
def save_results(is_running, _):
    if not is_running:
        raise dash.exceptions.PreventUpdate

    logs = []
    payload = cache.get(user_cache_key("user_result_array"))
    if not payload or not isinstance(payload, list) or len(payload) < 12:
        logs.append("❗ Cannot save results: dataset is empty or invalid.")
        return dash.no_update, logs, False, False

    import tempfile, os
    from datetime import datetime

    try:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".h5")
        tmp.close()
        # Ensure 14 slots + pad tags for older datasets
        if len(payload) < 14:
            payload = payload + [[] for _ in range(14 - len(payload))]
        ports_n = len(payload[1]) if isinstance(payload[1], list) else 0
        payload[12] = (payload[12] if len(payload) > 12 and isinstance(payload[12], list) else [])
        payload[13] = (payload[13] if len(payload) > 13 and isinstance(payload[13], list) else [])
        payload[12] = (payload[12] + [""] * ports_n)[:ports_n]
        payload[13] = (payload[13] + [""] * ports_n)[:ports_n]
        save_data_full_h5(tmp.name, payload)
        with open(tmp.name, "rb") as f:
            blob = f.read()
        os.remove(tmp.name)
        filename = f"results_{datetime.now().strftime('%Y-%m-%d_%H-%M')}.h5"
    except Exception as e:
        logs.append(f"❌ Failed to save results: {e}")
        return dash.no_update, logs, False, False

    logs.append(f"✅ Results saved as {filename}")
    return dcc.send_bytes(blob, filename=filename), logs, False, False


# ==========================================================
# 🟢 Python callback: Upload results (.h5) → cache (full 12-slot data_full)
# ==========================================================
@app.callback(
    Output("log-list-store", "data", allow_duplicate=True),
    Output("build-graph-btn", "disabled", allow_duplicate=True),
    Output("build-standard-btn", "disabled", allow_duplicate=True),
    Output("add-port-build-btn", "disabled", allow_duplicate=True),
    Output("precalc-spinner-output", "children", allow_duplicate=True),
    Output("build-standard-output", "children", allow_duplicate=True),
    Output("console-header", "children", allow_duplicate=True),
    Output("cache-update-trigger", "data", allow_duplicate=True),  # ensure allow_duplicate=True
    Input("upload-results", "contents"),
    State("upload-results", "filename"),
    State("cache-update-trigger", "data"),
    prevent_initial_call=True
)
def upload_results_to_cache(contents, filename, trgr):
    if not contents or not filename:
        raise dash.exceptions.PreventUpdate

    logs = []
    name = os.path.basename(filename)
    lower = name.lower()

    # Expecting HDF5
    if not (lower.endswith(".h5") or lower.endswith(".hdf5")):
        logs.append(f"❗ Unsupported file type: {name}. Please upload a .h5 file.")
        return logs, False, False, False, None, None, dash.no_update, trgr

    # Decode data URL

    try:
        header, b64 = contents.split(",", 1)
        if not header.startswith("data:"):
            raise ValueError("Invalid content header")
        decoded = base64.b64decode(b64, validate=True)
    except Exception as e:
        logs.append(f"❌ Failed to decode uploaded file: {e}")
        return logs, False, False, False, None, None, dash.no_update, trgr

    # Write to temp file then load full data_full
    try:
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".h5")
        tmp.write(decoded)
        tmp.close()
        data_full = load_data_full_h5(tmp.name)
        try:
            os.remove(tmp.name)
        except Exception:
            pass
    except Exception as e:
        logs.append(f"❌ Failed to parse HDF5: {e}")
        return logs, False, False, False, None, None, dash.no_update, trgr

    # Put into cache as full in-memory data_full (like Precalculate)
    cache.set(user_cache_key("user_result_array"), data_full)

    # Update header (uses your formatter internally)
    try:
        header = build_current_config_message(
            port_names=data_full[1],
            device_names=data_full[0],
            diag_types=data_full[11],
            config_tags=data_full[12],
            b0_tags=data_full[13],
        )
    except Exception:
        header = dash.no_update

    # Logs
    n_ports = len(data_full[1]) if isinstance(data_full[1], list) else 0
    logs.append(f"✅ Results loaded: {name}")
    logs.append(f"📦 Entries restored: {n_ports} ports")
    logs.append("🗄️ Cached in memory. You can Build/Add Port now.")

    return logs, False, False, False, None, None, header, (trgr or 0) + 1


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
