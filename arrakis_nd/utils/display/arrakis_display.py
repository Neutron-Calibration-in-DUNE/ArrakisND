import os
import yaml
from dash import (
    Dash,
    dcc,
    html,
    callback_context
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
        self.service_prefix = config['service_prefix']
        self.server_url = config['server_url']
        self.port = config['port']
        self.jupyter_mode = config['jupyter_mode']

        self.flow_folder = ''
        self.arrakis_folder = ''
        self.flow_files = []
        self.flow_file = ''
        self.arrakis_file = ''
        self.arrakis_files = []
        self.available_events = []
        self.unique_events = []
        self.event = None
        self.available_hits = []
        self.hit = None
        self.available_track_ids = []
        self.track_id = None
        self.hit_type = 'prompt'

        self.nersc_flow_folder = '/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/'
        self.standard_flow_folders = [
            {'label': 'MiniRun4', 'value':
                f'{self.nersc_flow_folder}MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW'},
            {'label': 'MiniRun4.5', 'value':
                f'{self.nersc_flow_folder}MiniRun4.5_1E19_RHC/inputs_beta3/MiniRun4.5_1E19_RHC.flow/FLOW/0000000'},
            {'label': 'MiniRun5', 'value':
                f'{self.nersc_flow_folder}MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.flow.beta2a/FLOW/0000000'},
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

        self.interactions = None
        self.segments = None
        self.stack = None
        self.trajectories = None
        self.charge = None
        self.light = None
        self.topology = None
        self.physics = None
        self.particle = None
        self.unique_topology = None

        self.left_tpc = TPCDisplay(id_suffix='left')
        self.right_tpc = TPCDisplay(id_suffix='right')
        self.left_larpix_light_display = LArPixLightDisplay()
        self.right_larpix_light_display = LArPixLightDisplay()

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
                # Right side: Menu items
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
        """Define the sidebar"""
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
                html.P("🔍 FLOW/ARRAKIS Folders & Files "),
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
                html.H2(),

                # Dropdown which lists available flow and arrakis files
                html.Label('FLOW files'),
                dcc.Dropdown(
                    id='flow_dropdown',
                    searchable=True,
                    placeholder='Select a FLOW file...',
                    style={'color': "#000000"}
                ),
                html.Label('ARRAKIS files'),
                dcc.Dropdown(
                    id='arrakis_dropdown',
                    searchable=True,
                    placeholder='Select an ARRAKIS file...',
                    style={'color': "#000000"}
                ),
                html.H2(),

                # Event selector and previous/next buttons
                html.H2(),
                html.Label('Event (spill)'),
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
                html.Div(id='event_output'),
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
            ],
            style=self.styles['SIDEBAR_STYLE'],
        )

    def construct_main_display(self):
        # Define content for the tabs
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
                                    {'label': 'TPC', 'value': 'tpc'},
                                    {'label': 'LArPix/Light', 'value': 'light'},
                                ],
                                value='tpc',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='left_window_plottype_dropdown',
                                options=[
                                    {'label': 'Q (charge)', 'value': 'q'},
                                    {'label': 'Topology', 'value': 'topology'},
                                    {'label': 'Physics', 'value': 'physics'},
                                    {'label': 'Particle', 'value': 'particle'},
                                    {'label': 'Tracklette', 'value': 'tracklette'},
                                    {'label': 'Track', 'value': 'track'},
                                    {'label': 'Fragment', 'value': 'fragment'},
                                    {'label': 'Shower', 'value': 'shower'},
                                    {'label': 'Blip', 'value': 'blip'},
                                ],
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
                        html.Label('Hit [hit_id]'),
                        html.Div([
                            dcc.Dropdown(
                                id='left_hit_dropdown',
                                options=self.available_hits,
                                searchable=True,
                                placeholder="Select an hit...",
                                style={'color': "#000000", 'width': '60%'}
                            ),
                            html.Button('Previous', id='left_previous_hit', style={'width': '20%'}),
                            html.Button('Next', id='left_next_hit', style={'width': '20%'}),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
                        html.H2(),
                        html.Label('TrackID [track_id]'),
                        html.Div([
                            dcc.Dropdown(
                                id='left_track_id_dropdown',
                                options=self.available_track_ids,
                                searchable=True,
                                placeholder="Select a track_id...",
                                style={'color': "#000000", 'width': '60%'}
                            ),
                            html.Button('Previous', id='left_previous_track_id', style={'width': '20%'}),
                            html.Button('Next', id='left_next_track_id', style={'width': '20%'}),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
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
                            ],
                            id='left_bottom_text',
                            style={'padding': '10px', 'margin-top': '10px'}
                        ),
                    ]),
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
                                    {'label': 'TPC', 'value': 'tpc'},
                                    {'label': 'LArPix/Light', 'value': 'light'},
                                ],
                                value='tpc',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='right_window_plottype_dropdown',
                                options=[
                                    {'label': 'Q (charge)', 'value': 'q'},
                                    {'label': 'Topology', 'value': 'topology'},
                                    {'label': 'Physics', 'value': 'physics'},
                                    {'label': 'Particle', 'value': 'particle'},
                                    {'label': 'Tracklette', 'value': 'tracklette'},
                                    {'label': 'Track', 'value': 'track'},
                                    {'label': 'Fragment', 'value': 'fragment'},
                                    {'label': 'Shower', 'value': 'shower'},
                                    {'label': 'Blip', 'value': 'blip'},
                                ],
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
                        html.H2(),
                        html.Label('Hit [hit_id]'),
                        html.Div([
                            dcc.Dropdown(
                                id='right_hit_dropdown',
                                options=self.available_hits,
                                searchable=True,
                                placeholder="Select an hit...",
                                style={'color': "#000000", 'width': '60%'}
                            ),
                            html.Button('Previous', id='right_previous_hit', style={'width': '20%'}),
                            html.Button('Next', id='right_next_hit', style={'width': '20%'}),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
                        html.H2(),
                        html.Label('TrackID [track_id]'),
                        html.Div([
                            dcc.Dropdown(
                                id='right_track_id_dropdown',
                                options=self.available_track_ids,
                                searchable=True,
                                placeholder="Select a track_id...",
                                style={'color': "#000000", 'width': '60%'}
                            ),
                            html.Button('Previous', id='right_previous_track_id', style={'width': '20%'}),
                            html.Button('Next', id='right_next_track_id', style={'width': '20%'}),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px'}),
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
                            ],
                            id='right_bottom_text',
                            style={'padding': '10px', 'margin-top': '10px'}
                        ),
                    ]),
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
        # Callbacks to update the content based on which dropdown is selected
        @self.app.callback(
            Output("dynamic_left_content", "children"),
            [Input("left_window_dropdown", "value")]
        )
        def render_left_content(value):
            if value == "light":
                return self.left_larpix_light_display.layout
            elif value == "tpc":
                return self.left_tpc.layout
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The value {value} was not recognised..."),
                ],
                className="p-3 bg-light rounded-3",
            )

        @self.app.callback(
            [Output('left_bottom_text', 'children'),
             Output('left_hit_dropdown', 'value'),
             Output('left_track_id_dropdown', 'value')],
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
            pdg_id = point_data['customdata'][3]
            track_id = point_data['customdata'][4]
            self.track_id = track_id

            return (
                html.Div([
                    html.P(f'x: {x:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'y: {y:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'z: {z:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'Q: {Q:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'E: {E:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'pdg_id: {pdg_id}', style={'margin': '0', 'padding': '0'})
                ]),
                self.hit,
                self.track_id
            )

        @self.app.callback(
            [Output('right_bottom_text', 'children'),
             Output('right_hit_dropdown', 'value'),
             Output('right_track_id_dropdown', 'value')],
            Input('tpc_plot_right', 'clickData'),
            [State('tpc_plot_right', 'figure')]
        )
        def display_right_click_data(clickData, tpc_plot):
            if clickData is None:
                raise PreventUpdate
            point_data = clickData['points'][0]

            if 'customdata' not in point_data:
                raise PreventUpdate

            hit_id = point_data['customdata'][0]
            self.hit = hit_id
            # self.right_tpc.highlight_point(hit_id)

            x = point_data['x']
            y = point_data['y']
            z = point_data['z']
            Q = point_data['customdata'][1]
            E = point_data['customdata'][2]
            pdg_id = point_data['customdata'][3]
            track_id = point_data['customdata'][4]
            self.track_id = track_id

            return (
                html.Div([
                    html.P(f'x: {x:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'y: {y:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'z: {z:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'Q: {Q:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'E: {E:.3f}', style={'margin': '0', 'padding': '0'}),
                    html.P(f'pdg_id: {pdg_id}', style={'margin': '0', 'padding': '0'})
                ]),
                self.hit,
                self.track_id
            )

        @self.app.callback(
            Output("dynamic_right_content", "children"),
            [Input("right_window_dropdown", "value")]
        )
        def render_right_content(value):
            if value == "light":
                return self.right_larpix_light_display.layout
            elif value == "tpc":
                return self.right_tpc.layout
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The value {value} was not recognised..."),
                ],
                className="p-3 bg-light rounded-3",
            )

        @self.app.callback(
            Output('flow_folder_input', 'value'),
            Input('standard_flow_dropdown', 'value')
        )
        def update_flow_folder(
            flow_folder
        ):
            return flow_folder

        # Callback to update dropdown options
        @self.app.callback(
            Output('flow_dropdown', 'options'),
            Input('flow_folder_input', 'value'),
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
                return flow_options
            return []

        # Callback to update dropdown options
        @self.app.callback(
            Output('arrakis_dropdown', 'options'),
            Input('arrakis_folder_input', 'value'),
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
                return arrakis_options
            return []

        @self.app.callback(
            Output('event_dropdown', 'options'),
            [Input('flow_dropdown', 'value'), Input('arrakis_dropdown', 'value')],
        )
        def update_available_events(flow_file, arrakis_file):
            self.available_events = []
            if flow_file is not None:
                try:
                    self.flow_file = flow_file
                    with h5py.File(self.flow_folder + flow_file, "r") as flow_file:
                        trajectories = flow_file['mc_truth/trajectories/data']
                        events = trajectories['event_id']
                        for key in self.geometry_info.keys():
                            try:
                                self.geometry_info[key] = flow_file[f'geometry_info/{key}/data'][:]
                            except Exception:
                                print(f"Issue with getting {key} from geometry_info")
                        self.left_tpc.set_geometry_info(self.geometry_info)
                        self.right_tpc.set_geometry_info(self.geometry_info)
                        self.left_larpix_light_display.set_geometry_info(self.geometry_info)
                        self.right_larpix_light_display.set_geometry_info(self.geometry_info)
                        self.unique_events = np.unique(events)
                        self.available_events = [
                            {'label': event, 'value': event}
                            for event in self.unique_events
                        ]
                except Exception:
                    pass
            if arrakis_file is not None:
                try:
                    self.arrakis_file = arrakis_file
                    with h5py.File(self.arrakis_folder + arrakis_file, "r") as arrakis_file:
                        pass
                except Exception:
                    pass
            self.event = None
            return self.available_events

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

        @self.app.callback(
            [Output('event_output', 'children'),
             Output('tpc_plot_left', 'figure'),
             Output('tpc_plot_right', 'figure'),
             Output('left_hit_dropdown', 'options'),
             Output('left_track_id_dropdown', 'options')],
            [Input('event_dropdown', 'value'),
             Input('left_window_plottype_dropdown', 'value'),
             Input('right_window_plottype_dropdown', 'value')],
        )
        def load_event(event, left_plottype, right_plottype):
            print_output = ''
            self.left_tpc.update_plottype(left_plottype)
            self.right_tpc.update_plottype(right_plottype)
            if event is not None:
                if self.flow_file:
                    with h5py.File(self.flow_folder + self.flow_file, "r") as flow_file:
                        interactions_events = flow_file['mc_truth/interactions/data']['event_id']
                        segments_events = flow_file['mc_truth/segments/data']['event_id']
                        stack_events = flow_file['mc_truth/stack/data']['event_id']
                        trajectories_events = flow_file['mc_truth/trajectories/data']['event_id']
                        charge_segments = flow_file[f'mc_truth/calib_{self.hit_type}_hit_backtrack/data']['segment_ids'].astype(int)
                        charge_fraction = flow_file[f'mc_truth/calib_{self.hit_type}_hit_backtrack/data']['fraction']
                        charge_fraction_mask = (charge_fraction == 0)
                        charge_segments[charge_fraction_mask] = -1
                        non_zero_charge_segments = [row[row != 0] for row in charge_fraction]
                        max_length = len(max(non_zero_charge_segments, key=len))
                        segments_ids = flow_file['mc_truth/segments/data']['segment_id']
                        self.interactions = flow_file['mc_truth/interactions/data'][
                            np.where(interactions_events == event)[0]
                        ]
                        self.segments = flow_file['mc_truth/segments/data'][
                            np.where(segments_events == event)[0]
                        ]
                        self.stack = flow_file['mc_truth/stack/data'][
                            np.where(stack_events == event)[0]
                        ]
                        self.trajectories = flow_file['mc_truth/trajectories/data'][
                            np.where(trajectories_events == event)[0]
                        ]
                        """For charge data we must backtrack through segments"""
                        hits_to_segments = np.any(
                            np.isin(
                                charge_segments[:, :max_length], segments_ids[(segments_events == event)]
                            ),
                            axis=1,
                        )
                        self.charge = flow_file[f'charge/calib_{self.hit_type}_hits/data'][
                            hits_to_segments
                        ]
                        self.available_hits = [ii for ii in range(len(self.charge))]
                        # charge_events = flow_file["charge/events/data"]["id"]

                        # self.charge_events = flow_file["charge/events/data"][np.where(interactions_events == event)[0]]

                        # """Likewise for light data, we must backtrack through segments"""
                        # match_light = flow_file['/light/events/data'][:][
                        #     flow_file['/charge/events/ref/light/events/ref'][np.where(interactions_events == event)[0], 1]
                        # ]["id"]
                        # we have to try them all, events may not be time ordered
                        # self.light = flow_file["light/events/data"][:]

                        # waveforms_all_detectors = flow_file["light/wvfm/data"]["samples"][match_light]
                        # print(match_light)
                        # we have now the waveforms for all detectors matched in time to the event
                    print_output = f"Event: {event} Loaded!"
                    self.left_tpc.update_flow_event(
                        self.interactions,
                        self.segments,
                        self.stack,
                        self.trajectories,
                        self.charge,
                    )
                    self.right_tpc.update_flow_event(
                        self.interactions,
                        self.segments,
                        self.stack,
                        self.trajectories,
                        self.charge,
                    )
                    # self.charge_light_display.construct_light_detectors(waveforms_all_detectors)
                    # self.charge_light_display.waveforms = self.charge_light_display.construct_waveforms()
                if self.arrakis_file:
                    with h5py.File(self.arrakis_folder + self.arrakis_file, "r") as arrakis_file:
                        arrakis_event_ids = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"]["event_id"]
                        self.topology = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"]["topology"][
                            (arrakis_event_ids == event)
                        ]
                        self.physics = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"]["physics"][
                            (arrakis_event_ids == event)
                        ]
                        self.particle = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"]["particle"][
                            (arrakis_event_ids == event)
                        ]
                        self.unique_topology = arrakis_file[f"charge/calib_{self.hit_type}_hits/data"]["unique_topology"][
                            (arrakis_event_ids == event)
                        ]
                        self.available_track_ids = np.unique(self.unique_topology)
                    self.left_tpc.update_arrakis_event(
                        self.topology,
                        self.physics,
                        self.particle,
                        self.unique_topology,
                    )
                    self.right_tpc.update_arrakis_event(
                        self.topology,
                        self.physics,
                        self.particle,
                        self.unique_topology
                    )
            self.left_tpc.plot_event()
            self.right_tpc.plot_event()
            return (
                print_output,
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
