import os
import yaml
from dash import (
    Dash,
    dcc,
    html,
    callback_context,
    no_update,
)
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import numpy as np
import h5py
import glob

from arrakis_nd.utils.display.set_server import SetServer
from arrakis_nd.utils.display.tpc_display import TPCDisplay
from arrakis_nd.utils.display.larpix_light_display import LArPixLightDisplay
from arrakis_nd.dataset.common import (
    process_type_dict,
    sub_process_type_dict
)


class ArrakisDisplay:
    """
    A class to create a Dash app for displaying the results of the ArrakisND analysis.
    """
    def __init__(
        self,
    ):
        """
        Initialize the ArrakisDisplay class.
        """

        server = SetServer()
        config = server.get_config()

        """Server parameters"""
        self.service_prefix = config['service_prefix']
        self.server_url = config['server_url']
        self.port = config['port']
        self.jupyter_mode = config['jupyter_mode']

        """Flow parameters"""
        self.flow_folder = ''
        self.flow_files = []
        self.flow_file = ''
        """Arrakis parameters"""
        self.arrakis_folder = ''
        self.arrakis_files = []
        self.arrakis_file = ''
        """Blip parameters"""
        self.blip_folder = ''
        self.blip_files = []
        self.blip_file = ''

        """Event information"""
        self.available_events = []
        self.unique_events = []
        self.event = None
        self.start_indices = []
        self.end_indices = []
        self.start_indices_map = {}
        self.end_indices_map = {}

        """Hits and TrackIDs"""
        self.available_hits = []
        self.hit = None
        self.available_track_ids = []
        self.track_id = None
        self.hit_type = 'prompt'

        """Standard NERSC Flow folders"""
        self.nersc_flow_folder = '/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/'
        self.standard_flow_folders = [
            {'label': 'MiniRun4', 'value':
                f'{self.nersc_flow_folder}MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW'},
            {'label': 'MiniRun4.5', 'value':
                f'{self.nersc_flow_folder}MiniRun4.5_1E19_RHC/inputs_beta3/MiniRun4.5_1E19_RHC.flow/FLOW/0000000'},
            {'label': 'MiniRun5', 'value':
                f'{self.nersc_flow_folder}MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.flow.beta2a/FLOW/0000000'},
            {'label': 'Half field data', 'value':
                ''},
        ]

        self.geometry_info = {
            'anode_drift_coordinate': [],
            'det_bounds': [],
            'det_id': [],
            'det_rel_pos': [],
            'drift_dir': [],
            'pixel_coordinates_2D': [],
            'sipm_abs_pos': [],
            'sipm_rel_pos': [],
            'tile_id': [],
        }

        """Data objects from flow/arrakis/blip"""
        """Flow"""
        self.flow_plot_options = [
            {'label': 'Q (charge)', 'value': 'q'},
            {'label': 'P (light)', 'value': 'p'},
        ]
        self.flow_truth = {
            'interactions': None,
            'segments': None,
            'stack': None,
            'trajectories': None
        }
        self.flow_event = {
            'charge': None,
            'light': None,
        }
        """Arrakis"""
        self.arrakis_plot_options = [
            {'label': 'Topology', 'value': 'topology'},
            {'label': 'Physics', 'value': 'physics'},
            {'label': 'Particle', 'value': 'particle'},
            {'label': 'Tracklette', 'value': 'tracklette'},
            {'label': 'Track', 'value': 'track'},
            {'label': 'Fragment', 'value': 'fragment'},
            {'label': 'Shower', 'value': 'shower'},
            {'label': 'Blip', 'value': 'blip'},
        ]
        self.arrakis_event = {
            'topology': None,
            'physics': None,
            'particle': None,
            'unique_topology': None,
            'vertex': None,
            'tracklette_begin': None,
            'tracklette_end': None,
            'fragment_begin': None,
            'fragment_end': None,
            'shower_begin': None,
        }
        """Blip"""
        self.blip_plot_options = [
            {'label': 'Topology', 'value': 'topology'},
            {'label': 'Physics', 'value': 'physics'},
            {'label': 'Particle', 'value': 'particle'},
            {'label': 'Tracklette', 'value': 'tracklette'},
            {'label': 'Track', 'value': 'track'},
            {'label': 'Fragment', 'value': 'fragment'},
            {'label': 'Shower', 'value': 'shower'},
            {'label': 'Blip', 'value': 'blip'},
        ]
        self.blip_event = {
            'topology': None,
            'physics': None,
            'particle': None,
            'unique_topology': None,
            'vertex': None,
            'tracklette_begin': None,
            'tracklette_end': None,
            'fragment_begin': None,
            'fragment_end': None,
            'shower_begin': None,
        }

        """TPC Displays"""
        self.left_tpc = TPCDisplay(id_suffix='left')
        self.right_tpc = TPCDisplay(id_suffix='right')
        """Light and LArPix Displays"""
        self.left_larpix_light_display = LArPixLightDisplay()
        self.right_larpix_light_display = LArPixLightDisplay()

        """Datapoint information"""
        self.empty_point_text = html.Div([
            html.P('x: ', style={'margin': '0', 'padding': '0'}),
            html.P('y: ', style={'margin': '0', 'padding': '0'}),
            html.P('z: ', style={'margin': '0', 'padding': '0'}),
            html.P('Q: ', style={'margin': '0', 'padding': '0'}),
            html.P('E: ', style={'margin': '0', 'padding': '0'}),
            html.P('pdg_id: ', style={'margin': '0', 'padding': '0'}),
        ])

        self.construct_app()
        self.construct_layout()
        self.construct_callbacks()
        self.run_app()

    def adjust_iframe_height(self, height=1000):
        """
        Generates a script to adjust the iframe height for the Dash app when running in Jupyter.
        Parameters:
        height (int): The desired height of the iframe in pixels.
        """
        from IPython.display import display, HTML
        script = f"""
        <script>
        // You might need to adjust the selector depending on your Jupyter environment
        const iframes = document.querySelectorAll('iframe');
        iframes.forEach(function(iframe) {{
            iframe.style.height = '{height}px';
        }});
        </script>
        """
        display(HTML(script))

    def construct_app(self):
        """
        Construct the Dash app.
        """
        if self.jupyter_mode == "inline":
            self.app = Dash(
                __name__,
                requests_pathname_prefix=f"{self.service_prefix}proxy/{self.port}/",
                external_stylesheets=[dbc.themes.SPACELAB]
            )
        else:
            self.app = Dash(
                __name__,
                external_stylesheets=[dbc.themes.SPACELAB]
            )

        """Get the custom style file"""
        with open('arrakis_nd/utils/display/assets/styles.yaml', 'r') as file:
            self.styles = yaml.safe_load(file)

    def construct_layout(self):
        """
        Construct the layout.
        """
        """Define the navbar with a dropdown"""
        self.construct_navbar()
        self.construct_sidebar()
        self.construct_main_display()

    def construct_navbar(self):
        """
        The navbar contains a clickable link to the ArrakisND repository,
        as well as links to various 2x2 sites.
        """
        self.navbar = html.Div(
            children=[
                html.A(
                    href="https://github.com/Neutron-Calibration-in-DUNE/ArrakisND",
                    target="_blank",  # Opens the link in a new tab
                    children=[
                        html.Img(
                            src='/assets/github-mark.png',
                            style={'height': '35px', 'marginRight': '15px'}
                        ),
                    ],
                    style={'display': 'flex', 'alignItems': 'center', 'textDecoration': 'none'},
                ),
                dbc.DropdownMenu(
                    [
                        dbc.DropdownMenuItem(
                            "2x2 sim (github)",
                            href="https://github.com/DUNE/2x2_sim",
                            external_link=True
                        ),
                        dbc.DropdownMenuItem(
                            "MiniRun5 File Locations",
                            href="https://github.com/DUNE/2x2_sim/wiki/MiniRun5-file-locations",
                            external_link=True
                        ),
                    ],
                    label="Menu",
                ),
            ],
            style=self.styles['NAVBAR_STYLE']
        )

    def construct_sidebar(self):
        """
        The sidebar contains a set of text boxes for inputing
        directory and file information for FLOW/ARRAKIS/BLIP.
        Once a flow or arrakis or blip file is selected, the
        available events from that file are populated.
        """
        self.sidebar = html.Div(
            [
                # DUNE logo and Arrakis Display text
                html.Div(
                    children=[
                        html.Img(
                            src='/assets/2x2.png',
                            style={'height': '100px', 'marginRight': '15px'}
                        ),
                        html.H3("2x2 Display"),
                    ],
                    style={'display': 'flex', 'alignItems': 'center'}
                ),

                # Standard FLOW folder selection dropdown
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                html.P("🔍 FLOW Folders & Files "),
                dcc.Dropdown(
                    id='standard_flow_dropdown',
                    options=self.standard_flow_folders,
                    style={'color': "#000000"},
                ),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                # Text box for writing FLOW and Arrakis folders
                dbc.Label("🗃️ FLOW folder"),
                dbc.Input(
                    placeholder="Enter the FLOW folder",
                    type="text",
                    id='flow_folder_input',
                    size='sm',
                    value=''
                ),
                dbc.Label("🗃️ ARRAKIS folder"),
                dbc.Input(
                    placeholder="Enter the ARRAKIS folder",
                    type="text",
                    id='arrakis_folder_input',
                    size='sm',
                    value=''
                ),
                dbc.Label("🗃️ BLIP folder"),
                dbc.Input(
                    placeholder="Enter the BLIP folder",
                    type="text",
                    id='blip_folder_input',
                    size='sm',
                    value=''
                ),
                html.H2(),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                # Dropdown which lists available flow and arrakis files
                html.Label('FLOW files:'),
                dcc.Dropdown(
                    id='flow_dropdown',
                    searchable=True,
                    placeholder='Select a FLOW file...',
                    style={'color': "#000000"}
                ),
                html.Label('ARRAKIS files:'),
                dcc.Dropdown(
                    id='arrakis_dropdown',
                    searchable=True,
                    placeholder='Select an ARRAKIS file...',
                    style={'color': "#000000"}
                ),
                html.Label('BLIP files:'),
                dcc.Dropdown(
                    id='blip_dropdown',
                    searchable=True,
                    placeholder='Select an BLIP file...',
                    style={'color': "#000000"}
                ),
                html.H2(),
                # Event selector and previous/next buttons
                html.H2(),
                html.Label('Event (spill):'),
                html.Div([
                    dcc.Dropdown(
                        id='event_dropdown',
                        options=self.available_events,
                        searchable=True,
                        placeholder="Select an event...",
                        style={'color': "#000000", 'width': '60%'}
                    ),
                    html.Button('Previous', id='previous_event', style={'width': '20%'}),
                    html.Button('Next', id='next_event', style={'width': '20%'}),
                ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
                html.H2(),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                html.Label('Display status:'),
                html.Div(id='display_status', style={
                    'whiteSpace': 'normal',  # Allow text to wrap),
                }),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                html.H2(),
                html.Label("Hit type"),
                dcc.Dropdown(
                    id='hit_type_dropdown',
                    options=[
                        {'label': 'prompt', 'value': 'prompt'},
                        {'label': 'final', 'value': 'final'},
                    ],
                    value='prompt',
                    style={'color': "#000000"}
                ),
                html.H2(),
                html.Label('Hit [hit_id]'),
                html.Div([
                    dcc.Dropdown(
                        id='hit_dropdown',
                        options=self.available_hits,
                        searchable=True,
                        placeholder="Select an hit...",
                        style={'color': "#000000", 'width': '60%'}
                    ),
                    html.Button('Previous', id='previous_hit', style={'width': '20%'}),
                    html.Button('Next', id='next_hit', style={'width': '20%'}),
                ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
                html.H2(),
                html.Label('TrackID [track_id]'),
                html.Div([
                    dcc.Dropdown(
                        id='track_id_dropdown',
                        options=self.available_track_ids,
                        searchable=True,
                        placeholder="Select a track_id...",
                        style={'color': "#000000", 'width': '60%'}
                    ),
                    html.Button('Previous', id='previous_track_id', style={'width': '20%'}),
                    html.Button('Next', id='next_track_id', style={'width': '20%'}),
                ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
            ],
            style=self.styles['SIDEBAR_STYLE'],
        )

    def construct_main_display(self):
        """
        The main display contains two windows which can be selected as
        different plot types.
        """
        self.main_display = html.Div(id="page_content", style=self.styles["CONTENT_STYLE"], children=[
            html.Div(id="left_window", style=self.styles["LEFT_COLUMN"], children=[
                html.Div(id="dynamic_left_content")
            ]),
            html.Div(id="right_window", style=self.styles["RIGHT_COLUMN"], children=[
                html.Div(id="dynamic_right_content")
            ]),
            html.Div(id="bottom_section", style=self.styles["BOTTOM_SECTION"], children=[
                html.Div(id="left_bottom_section", style=self.styles["LEFT_BOTTOM_SECTION"], children=[
                    html.Div([
                        html.Div(style={'display': 'flex', 'width': '100%', 'gap': '10px'}, children=[
                            html.Div(children=[
                                html.H2(),
                                html.Label('Window I Display type')
                            ], style={'width': '50%'}),
                            html.Div(children=[
                                html.H2(),
                                html.Label('Plot type (labels)')
                            ], style={'width': '50%'}),
                        ]),
                        html.Div([
                            dcc.Dropdown(
                                id='left_window_dropdown',
                                options=[
                                    {'label': 'FLOW', 'value': 'flow'},
                                    {'label': 'ARRAKIS', 'value': 'arrakis'},
                                    {'label': 'BLIP', 'value': 'blip'}
                                ],
                                value='flow',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='left_window_plottype_dropdown',
                                options=self.flow_plot_options,
                                value='q',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                        ], style={
                            'display': 'flex',
                            'flexDirection': 'row',
                            'gap': '10px', 'width': '100%'
                        }),
                        html.H2(),
                        html.Label('Segment/Hit Scale'),
                        dcc.Slider(
                            id='left_window_scale_slider',
                            min=0.0,
                            max=1.0,
                            step=0.01,
                            value=0.01,  # Default scale
                            marks={i / 10.0: f'{i / 10.0}' for i in range(0, 11)},
                        ),
                        html.H2(),
                        html.Label('Show Keypoints'),

                    ], style={'width': '50%'}),
                    html.Div([
                        html.H2(),
                        html.Label('Segment/Hit Info'),
                        html.Div(
                            [
                                html.P('x: ', style={'margin': '0', 'padding': '0'}),
                                html.P('y: ', style={'margin': '0', 'padding': '0'}),
                                html.P('z: ', style={'margin': '0', 'padding': '0'}),
                                html.P('Q: ', style={'margin': '0', 'padding': '0'}),
                                html.P('E: ', style={'margin': '0', 'padding': '0'}),
                                html.P('pdg_id: ', style={'margin': '0', 'padding': '0'}),
                                html.P('start_process: ', style={'margin': '0', 'padding': '0'}),
                                html.P('start_subprocess: ', style={'margin': '0', 'padding': '0'}),
                                html.P('parent_track_id: ', style={'margin': '0', 'padding': '0'}),
                                html.P('parent_pdg_id: ', style={'margin': '0', 'padding': '0'}),
                                html.P('parent_start_process: ', style={'margin': '0', 'padding': '0'}),
                                html.P('parent_start_subprocess: ', style={'margin': '0', 'padding': '0'}),
                            ],
                            id='bottom_text',
                            style={'padding': '10px', 'margin-top': '10px'}
                        ),
                    ], style={'width': '50%'}),
                ]),
                html.Div(id="right_bottom_section", style=self.styles["RIGHT_BOTTOM_SECTION"], children=[
                    html.Div([
                        html.Div(style={'display': 'flex', 'width': '100%', 'gap': '10px'}, children=[
                            html.Div(children=[
                                html.H2(),
                                html.Label('Window II Display type')
                            ], style={'width': '50%'}),
                            html.Div(children=[
                                html.H2(),
                                html.Label('Plot type (labels)')
                            ], style={'width': '50%'}),
                        ]),
                        html.Div([
                            dcc.Dropdown(
                                id='right_window_dropdown',
                                options=[
                                    {'label': 'FLOW', 'value': 'flow'},
                                    {'label': 'ARRAKIS', 'value': 'arrakis'},
                                    {'label': 'BLIP', 'value': 'blip'}
                                ],
                                value='flow',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='right_window_plottype_dropdown',
                                options=self.flow_plot_options,
                                value='q',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                        ], style={
                            'display': 'flex',
                            'flexDirection': 'row',
                            'gap': '10px', 'width': '100%'
                        }),
                        html.H2(),
                        html.Label('Segment/Hit Scale'),
                        dcc.Slider(
                            id='right_window_scale_slider',
                            min=0.0,
                            max=1.0,
                            step=0.01,
                            value=0.01,  # Default scale
                            marks={i / 10.0: f'{i / 10.0}' for i in range(0, 11)},
                        ),
                    ], style={'width': '50%'}),
                ]),
            ])
        ])

        # Define the layout
        self.app.layout = html.Div(style={'overflow': 'scroll'}, children=[
            dcc.Location(id="url"),
            self.navbar,
            self.sidebar,
            self.main_display,
            dcc.Store(id='no_output', data=0)
        ])

    def construct_callbacks(self):
        """
        Standard FLOW folder dropdown.  This gives a few options that
        can be selected to automatically fill the flow folder text input
        with different MiniRun locations.
        """
        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('flow_folder_input', 'value')],
            Input('standard_flow_dropdown', 'value'),
            prevent_initial_call=True
        )
        def update_flow_folder(
            flow_folder
        ):
            if flow_folder is not None:
                print_status = 'Setting flow folder'
            else:
                print_status = ''
            return (print_status, flow_folder)

        """
        FLOW, ARRAKIS and BLIP folder and file inputs. These callbacks
        take changes from the flow/arrakis/blip folder inputs and
        automatically search recursively for the respective hdf5 files
        within those folders.  The file callback loads the files and
        determines what events are present.  Those events are then populated
        in the event selector.
        """
        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('flow_dropdown', 'options')],
            Input('flow_folder_input', 'value'),
            prevent_initial_call=True
        )
        def update_flow_folder_files(
            flow_folder
        ):
            """Check that flow folder has a '/' at the end"""
            if flow_folder:
                if flow_folder[-1] != '/':
                    flow_folder += '/'
            self.flow_folder = flow_folder

            flow_options = []
            if flow_folder and os.path.isdir(flow_folder):
                self.flow_files = sorted([
                    os.path.basename(input_file) for input_file in glob.glob(
                        f"{flow_folder}*.hdf5", recursive=True
                    )
                    if 'FLOW' in input_file
                ])
                flow_options = [
                    {'label': file, 'value': file}
                    for file in self.flow_files
                ]
                return (
                    f'Found {len(flow_options)} FLOW files',
                    flow_options
                )
            return (
                '',
                []
            )

        # Callback to update dropdown options
        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('arrakis_dropdown', 'options')],
            Input('arrakis_folder_input', 'value'),
            prevent_initial_call=True,
        )
        def update_arrakis_folder_files(
            arrakis_folder
        ):
            """Check that arrakis folder has a '/' at the end"""
            if arrakis_folder:
                if arrakis_folder[-1] != '/':
                    arrakis_folder += '/'

            self.arrakis_folder = arrakis_folder

            arrakis_options = []
            if arrakis_folder and os.path.isdir(arrakis_folder):
                self.arrakis_files = sorted([
                    os.path.basename(input_file) for input_file in glob.glob(
                        f"{arrakis_folder}*.hdf5", recursive=True
                    )
                    if 'ARRAKIS' in input_file
                ])
                arrakis_options = [
                    {'label': file, 'value': file}
                    for file in self.arrakis_files
                ]
                return (
                    f'Found {len(arrakis_options)} ARRAKIS files',
                    arrakis_options
                )
            return (
                '',
                []
            )

        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('blip_dropdown', 'options')],
            Input('blip_folder_input', 'value'),
            prevent_initial_call=True,
        )
        def update_blip_folder_files(
            blip_folder
        ):
            """Check that blip folder has a '/' at the end"""
            if blip_folder:
                if blip_folder[-1] != '/':
                    blip_folder += '/'

            self.blip_folder = blip_folder

            blip_options = []
            if blip_folder and os.path.isdir(blip_folder):
                self.blip_files = sorted([
                    os.path.basename(input_file) for input_file in glob.glob(
                        f"{blip_folder}*.hdf5", recursive=True
                    )
                    if 'BLIP' in input_file
                ])
                blip_options = [
                    {'label': file, 'value': file}
                    for file in self.blip_files
                ]
                return(
                    f'Found {len(blip_options)} BLIP files',
                    blip_options
                )
            return (
                '',
                []
            )

        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('event_dropdown', 'options')],
            Input('flow_dropdown', 'value'),
            prevent_initial_call=True,
        )
        def update_available_events(flow_file):
            self.available_events = []
            display_status = ''
            if flow_file is not None:
                try:
                    self.flow_file = flow_file
                    with h5py.File(self.flow_folder + flow_file, "r") as flow_file:
                        try:
                            """Update the TPC objects with the interactions, segments, stacks and trajectories"""
                            for key in self.flow_truth.keys():
                                self.flow_truth[key] = flow_file[f'mc_truth/{key}/data'][:]
                        except Exception:
                            display_status += 'No truth info in FLOW file.'
                            for key in self.flow_truth.keys():
                                self.flow_truth[key] = None
                        """Send truth info to TPCs"""
                        self.left_tpc.update_flow_truth(
                            self.flow_truth
                        )
                        self.right_tpc.update_flow_truth(
                            self.flow_truth
                        )
                        try:
                            """Get event information from FLOW"""                            
                            events = flow_file['charge/events/data']
                            event_id = events['id']
                            nhits = events['nhit']
                            event_ids = []
                            for jj in range(len(event_id)):
                                event_ids += [event_id[jj] for kk in range(nhits[jj])]
                            self.unique_events, start_indices = np.unique(event_ids, return_index=True)
                            self.start_indices = start_indices.tolist()
                            self.end_indices = self.start_indices[1:] + [len(event_ids)]
                            self.start_indices_map = {
                                event: self.start_indices[jj]
                                for jj, event in enumerate(self.unique_events)
                            }
                            self.end_indices_map = {
                                event: self.end_indices[jj]
                                for jj, event in enumerate(self.unique_events)
                            }
                            display_status += f'Found {len(self.unique_events)} events in FLOW file.'
                        except Exception:
                            display_status += 'Issue getting event indices from flow file.'
                        for key in self.geometry_info.keys():
                            try:
                                self.geometry_info[key] = flow_file[f'geometry_info/{key}/data'][:]
                            except Exception:
                                display_status += f'Issue with getting {key} from geometry_info.'
                        self.left_tpc.set_geometry_info(self.geometry_info)
                        self.right_tpc.set_geometry_info(self.geometry_info)
                        self.left_larpix_light_display.set_geometry_info(self.geometry_info)
                        self.right_larpix_light_display.set_geometry_info(self.geometry_info)
                        self.available_events = [
                            {'label': event, 'value': event}
                            for event in self.unique_events
                        ]
                except Exception:
                    display_status += 'Issue loading FLOW file.'
            self.event = None
            return display_status, self.available_events

        @self.app.callback(
            Output('display_status', 'children', allow_duplicate=True),
            Input('arrakis_dropdown', 'value'),
            prevent_initial_call=True,
        )
        def update_arrakis_file(arrakis_file):
            display_status = ''
            if arrakis_file is not None:
                try:
                    self.arrakis_file = arrakis_file
                except Exception:
                    display_status = 'Issue setting arrakis file'
            return display_status
        
        @self.app.callback(
            Output('display_status', 'children', allow_duplicate=True),
            Input('blip_dropdown', 'value'),
            prevent_initial_call=True,
        )
        def update_blip_file(blip_file):
            display_status = ''
            if blip_file is not None:
                try:
                    self.blip_file = blip_file
                except Exception:
                    display_status = 'Issue setting blip file'
            return display_status

        @self.app.callback(
            Output('event_dropdown', 'value'),
            [
                Input('previous_event', 'n_clicks'),
                Input('next_event', 'n_clicks')
            ],
            [State('event_dropdown', 'value')]
        )
        def update_event(previous_clicks, next_clicks, current_value):
            triggered_id = callback_context.triggered[0]['prop_id'].split('.')[0] if callback_context.triggered else ''

            if not self.available_events:
                raise PreventUpdate

            current_index = next(
                (i for i, event in enumerate(self.available_events)
                 if event['value'] == current_value), None
            )

            if current_index is None:
                raise PreventUpdate

            new_index = current_index
            if triggered_id == 'previous_event' and current_index > 0:
                new_index = current_index - 1
            elif triggered_id == 'next_event' and current_index < len(self.available_events) - 1:
                new_index = current_index + 1
            else:
                # If we can't go previous or next, raise PreventUpdate to do nothing
                raise PreventUpdate
            self.event = self.available_events[new_index]['value']
            return self.available_events[new_index]['value']

        """
        Left and right window callbacks.  These are associated to the
        dropdowns that select the type of plot to show in the window.
        """
        @self.app.callback(
            [Output("dynamic_left_content", "children"),
             Output("left_window_plottype_dropdown", "options")],
            [Input("left_window_dropdown", "value")]
        )
        def render_left_content(value):
            self.left_tpc.update_datatype(value)
            if value == "flow":
                return self.left_tpc.layout, self.flow_plot_options
            elif value == "arrakis":
                return self.left_tpc.layout, self.arrakis_plot_options
            elif value == "blip":
                return self.left_tpc.layout, self.blip_plot_options
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The value {value} was not recognised..."),
                ],
                className="p-3 bg-light rounded-3",
            )

        @self.app.callback(
            [Output("dynamic_right_content", "children"),
             Output("right_window_plottype_dropdown", "options")],
            [Input("right_window_dropdown", "value")]
        )
        def render_right_content(value):
            self.right_tpc.update_datatype(value)
            if value == "flow":
                return self.right_tpc.layout, self.flow_plot_options
            elif value == "arrakis":
                return self.right_tpc.layout, self.arrakis_plot_options
            elif value == "blip":
                return self.right_tpc.layout, self.blip_plot_options
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The value {value} was not recognised..."),
                ],
                className="p-3 bg-light rounded-3",
            )

        """
        Left and right bottom text and hit/track_id dropdowns. These
        are associated to the bottom text information and the hit id
        and track id of selected points within a plot.
        """
        @self.app.callback(
            [Output('bottom_text', 'children'),
             Output('hit_dropdown', 'value'),
             Output('track_id_dropdown', 'value')],
            Input('tpc_plot_left', 'clickData'),
            [State('tpc_plot_left', 'figure')]
        )
        def display_left_click_data(clickData, tpc_plot):
            if clickData is None:
                raise PreventUpdate
            point_data = clickData['points'][0]

            if 'customdata' not in point_data:
                raise PreventUpdate

            hit_id = point_data['customdata'][0]
            self.hit = hit_id
            # self.left_tpc.highlight_point(hit_id)

            x = point_data['x']
            y = point_data['y']
            z = point_data['z']
            Q = point_data['customdata'][1]
            E = point_data['customdata'][2]
            try:
                pdg_id = point_data['customdata'][3]
                track_id = point_data['customdata'][4]
                self.track_id = track_id
                vertex_id = self.flow_truth['trajectories']['vertex_id'][track_id]
                parent_id = self.flow_truth['trajectories']['parent_id'][track_id]
                parent_index = np.where(
                    (self.flow_truth['trajectories']['traj_id'] == parent_id) & (self.flow_truth['trajectories']['vertex_id'] == vertex_id)
                )[0][0]
                start_process = process_type_dict[self.flow_truth['trajectories']['start_process'][track_id]]
                start_subprocess = sub_process_type_dict[self.flow_truth['trajectories']['start_subprocess'][track_id]]
                parent_pdg_id = self.flow_truth['trajectories']['pdg_id'][parent_index]
                parent_start_process = process_type_dict[self.flow_truth['trajectories']['start_process'][parent_index]]
                parent_start_subprocess = sub_process_type_dict[self.flow_truth['trajectories']['start_subprocess'][parent_index]]
            except Exception:
                pdg_id = -1
                track_id = -1
                self.track_id = -1
                vertex_id = -1
                parent_id = -1
                start_process = -1
                start_subprocess = -1
                parent_pdg_id = -1
                parent_start_process = -1
                parent_start_subprocess = -1

            return (
                html.Div([
                    html.P(f'x: {x:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'y: {y:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'z: {z:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'Q: {Q:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'E: {E:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'pdg_id: {pdg_id}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'start_process: {start_process}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'start_subprocess: {start_subprocess}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'parent_track_id: {parent_index}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'parent_pdg_id: {parent_pdg_id}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'parent_start_process: {parent_start_process}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'parent_start_subprocess: {parent_start_subprocess}', style={'margin': '0', 'padding': '0'}),
                ]),
                self.hit,
                self.track_id
            )

        """
        This callback updates the data for each plot type when a new event is loaded.
        The data is extracted from the flow/arrakis/blip files and sent to the different
        TPC Displays.
        """
        @self.app.callback(
            [Output('display_status', 'children', allow_duplicate=True),
             Output('tpc_plot_left', 'figure'),
             Output('tpc_plot_right', 'figure'),
             Output('hit_dropdown', 'options'),
             Output('track_id_dropdown', 'options')],
            [Input('event_dropdown', 'value'),
             Input('left_window_plottype_dropdown', 'value'),
             Input('right_window_plottype_dropdown', 'value')],
            prevent_initial_call=True
        )
        def load_event(event, left_plottype, right_plottype):
            display_output = ''
            try:
                self.left_tpc.update_plottype(left_plottype)
                self.right_tpc.update_plottype(right_plottype)
            except Exception as exception:
                display_output += f'Event loading ERROR: {exception}'
            if event is not None:
                display_output += f"Event: {event} Loaded!"
                try:
                    if self.flow_file:
                        with h5py.File(self.flow_folder + self.flow_file, "r") as flow_file:
                            self.flow_event['charge'] = flow_file[f'charge/calib_{self.hit_type}_hits/data'][
                                self.start_indices_map[event]:self.end_indices_map[event]
                            ]
                            self.available_hits = [ii for ii in range(len(self.flow_event['charge']))]
                        self.left_tpc.update_flow_event(
                            self.flow_event,
                        )
                        self.right_tpc.update_flow_event(
                            self.flow_event,
                        )
                        # self.charge_light_display.construct_light_detectors(waveforms_all_detectors)
                        # self.charge_light_display.waveforms = self.charge_light_display.construct_waveforms()
                except Exception as exception:
                    display_output += f'ERROR getting flow/arrakis: {exception}'
                try:
                    if self.arrakis_file:
                        with h5py.File(self.arrakis_folder + self.arrakis_file, "r") as arrakis_file:
                            for key in self.arrakis_event.keys():
                                try:
                                    self.arrakis_event[key] = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"][key][
                                        self.start_indices_map[event]:self.end_indices_map[event]
                                    ]
                                except Exception:
                                    display_output += f'ERROR getting {key} from ARRAKIS event.'
                                    self.arrakis_event[key] = None
                            self.available_track_ids = np.unique(self.arrakis_event['unique_topology'])
                        self.left_tpc.update_arrakis_event(
                            self.arrakis_event
                        )
                        self.right_tpc.update_arrakis_event(
                            self.arrakis_event
                        )
                except Exception as exception:
                    display_output += f'ERROR updating tpcs: {exception}'
                try:
                    if self.blip_file:
                        with h5py.File(self.blip_folder + self.blip_file, "r") as blip_file:
                            for key in self.blip_event.keys():
                                try:
                                    self.blip_event[key] = blip_file[f"charge/calib_{self.hit_type}_hits/data"][key][
                                        self.start_indices_map[event]:self.end_indices_map[event]
                                    ]
                                except Exception:
                                    display_output += f'ERROR getting {key} from BLIP event'
                                    self.blip_event[key] = None                        
                        self.left_tpc.update_blip_event(
                            self.blip_event
                        )
                        self.right_tpc.update_blip_event(
                            self.blip_event
                        ) 
                except Exception as exception:
                    display_output += f'BLIP loading ERROR: {exception}'
            try:
                self.left_tpc.plot_event()
            except Exception as exception:
                display_output += f'Left TPC plotting ERROR: {exception}'
            try:
                self.right_tpc.plot_event()
            except Exception as exception:
                display_output += f'Right TPC plotting ERROR: {exception}'
            return (
                display_output,
                self.left_tpc.tpc,
                self.right_tpc.tpc,
                self.available_hits,
                self.available_track_ids
            )

        @self.app.callback(
            Output('light-waveform', 'figure'),
            [Input('charge_light_tpc_plot', 'figure'),
             Input('event_dropdown', 'value'),
             Input('charge_light_tpc_plot', 'clickData')]
        )
        def update_light_waveform(charge_light_tpc_plot, event, click_data):
            if click_data:
                try:
                    if self.flow_file:
                        with h5py.File(self.flow_folder + self.flow_file, "r") as flow_file:
                            interactions_events = flow_file['mc_truth/interactions/data']['event_id']

                            """Likewise for light data, we must backtrack through segments"""
                            match_light = flow_file['/light/events/data'][:][
                                flow_file['/charge/events/ref/light/events/ref'][np.where(interactions_events == event)[0], 1]
                            ]["id"]

                            waveforms_all_detectors = flow_file["light/wvfm/data"]["samples"][match_light]
                            opid = click_data['points'][0]['id'].split('_')[1]
                            self.charge_light_display.plot_waveform(opid, waveforms_all_detectors)

                            print("waveform plotted")
                            print(opid)
                            print(self.charge_light_display.waveforms.data)

                            self.charge_light_display.generate_layout()
                    return self.charge_light_display.waveforms
                except Exception as e:
                    print(e)
                    print("that is not a light trap, no waveform to plot")

        @self.app.callback(
            Output('no_output', 'data'),
            Input('left_window_scale_slider', 'value'),
            [State('tpc_left_plot', 'figure')]
        )
        def update_left_size_scaler(size, tpc_plot):
            self.left_tpc.scale = size
            if self.event is not None:
                self.left_tpc.plot_event()
            return size

        # @self.app.callback(
        #     Output('tpc_plot_right', 'figure'),
        #     Input('right_window_scale_slider', 'value')
        # )
        # def update_right_size_scaler(size):
        #     print(size)
        #     self.right_tpc.scale = size
        #     if self.event is not None:
        #         self.right_tpc.plot_event()
        #     return self.right_tpc.tpc

    def run_app(self):
        self.app.run_server(
            jupyter_mode=self.jupyter_mode,
            jupyter_server_url=self.server_url,
            host="localhost",
            port=self.port,
        )
        if self.jupyter_mode == "inline":
            self.adjust_iframe_height(height=1500)


if __name__ == "__main__":
    display = ArrakisDisplay()
