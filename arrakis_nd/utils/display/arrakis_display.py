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

        self.left_tpc = TPCDisplay(id_suffix='left')
        self.right_tpc = TPCDisplay(id_suffix='right')
        self.left_larpix_light_display = LArPixLightDisplay()
        self.right_larpix_light_display = LArPixLightDisplay()

        self.construct_app()
        self.construct_widgets()
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
            styles = yaml.safe_load(file)

        """Define the navbar with a dropdown"""
        self.navbar = html.Div(
            children=[
                html.A(
                    href="https://github.com/Neutron-Calibration-in-DUNE/ArrakisND",
                    target="_blank",  # Opens the link in a new tab
                    children=[
                        html.Img(src='/assets/github-mark.png', style={'height': '35px', 'marginRight': '15px'}),
                    ],
                    style={'display': 'flex', 'alignItems': 'center', 'textDecoration': 'none'},
                ),
                # Right side: Menu items
                dbc.DropdownMenu(
                    [
                        dbc.DropdownMenuItem(
                            "A button", id="dropdown-button", n_clicks=0
                        ),
                        dbc.DropdownMenuItem(
                            "Internal link", href="/docs/components/dropdown_menu"
                        ),
                        dbc.DropdownMenuItem(
                            "External Link", href="https://github.com"
                        ),
                        dbc.DropdownMenuItem(
                            "External relative",
                            href="/docs/components/dropdown_menu",
                            external_link=True,
                        ),
                    ],
                    label="Menu",
                ),
            ],
            style=styles['NAVBAR_STYLE']
        )

        """Define the sidebar"""
        self.sidebar = html.Div(
            [
                # DUNE logo and Arrakis Display text
                html.Div(
                    children=[
                        html.Img(src='/assets/2x2.png', style={'height': '100px', 'marginRight': '15px'}),
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
                dbc.Label("FLOW folder"),
                dbc.Input(
                    placeholder="Enter the FLOW folder",
                    type="text",
                    id='flow_folder_input',
                    size='sm',
                    value=''
                ),
                dbc.Label("ARRAKIS folder"),
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

                # Highlighting options
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                # html.P("FLOW charge types to plot"),
                # dcc.Checklist(
                #     ['Segments', 'Calib Prompt Hits', 'Calib Final Hits'],
                #     ['Calib Final Hits'], inline=True
                # ),
                # # Scale Slider
                # html.Hr(style={'border': '1px solid #ffffff', 'height': '0px'}),
                # html.Div([
                #     html.Div([
                #         html.P("Active TPCs"),
                #         dcc.Checklist(
                #             ['0', '1', '2', '3', '4', '5', '6', '7'],
                #             ['0', '1', '2', '3', '4', '5', '6', '7']
                #         )], style={'width': '45vw'}),
                #     html.Div([
                #         html.P("Detectors"),
                #         dcc.Checklist(
                #             ['0', '1', '2', '3', '4', '5', '6', '7'],
                #             ['0', '1', '2', '3', '4', '5', '6', '7']
                #         )], style={'width': '45vw'}),
                # ], style={'display': 'flex'}),
            ],
            style=styles['SIDEBAR_STYLE'],
        )

        # Define content for the tabs
        self.main_display = html.Div(id="page-content", style=styles["CONTENT_STYLE"], children=[
            html.Div(id="left-window", style=styles["LEFT_COLUMN"], children=[
                html.Div(id="dynamic-left-content")
            ]),
            html.Div(id="right-window", style=styles["RIGHT_COLUMN"], children=[
                html.Div(id="dynamic-right-content")
            ]),
            html.Div(id="bottom-section", style=styles["BOTTOM_SECTION"], children=[
                html.Div(id="left-bottom-section", style=styles["LEFT_BOTTOM_SECTION"], children=[
                    html.Div([
                        html.Div(style={'display': 'flex', 'width': '100%', 'gap': '10px'}, children=[
                            html.Div(children=[html.H2(), html.Label('Window I Display type')], style={'width': '50%'}),
                            html.Div(children=[html.H2(), html.Label('Plot type (labels)')], style={'width': '50%'}),
                        ]),
                        html.Div([
                            dcc.Dropdown(
                                id='left-window-dropdown',
                                options=[
                                    {'label': 'TPC', 'value': 'tpc'},
                                    {'label': 'LArPix/Light', 'value': 'light'},
                                ],
                                value='tpc',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='left-window-plottype-dropdown',
                                options=[
                                    {'label': 'Q (charge)', 'value': 'q'},
                                    {'label': 'Topology', 'value': 'topology'},
                                    {'label': 'Physics', 'value': 'physics'},
                                ],
                                value='q',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px', 'width': '100%'}),
                        html.H2(),
                        html.Label('Segment/Hit Scale'),
                        dcc.Slider(
                            id='left-window-scale-slider',
                            min=0.0,
                            max=1.0,
                            step=0.01,
                            value=0.01,  # Default scale
                            marks={i / 10.0: f'{i / 10.0}' for i in range(0, 11)},
                        ),
                    ], style={'width': '50%'}),
                ]),
                html.Div(id="right-bottom-section", style=styles["RIGHT_BOTTOM_SECTION"], children=[
                    html.Div([
                        html.Div(style={'display': 'flex', 'width': '100%', 'gap': '10px'}, children=[
                            html.Div(children=[html.H2(), html.Label('Window II Display type')], style={'width': '50%'}),
                            html.Div(children=[html.H2(), html.Label('Plot type (labels)')], style={'width': '50%'}),
                        ]),
                        html.Div([
                            dcc.Dropdown(
                                id='right-window-dropdown',
                                options=[
                                    {'label': 'TPC', 'value': 'tpc'},
                                    {'label': 'LArPix/Light', 'value': 'light'},
                                ],
                                value='tpc',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                            dcc.Dropdown(
                                id='right-window-plottype-dropdown',
                                options=[
                                    {'label': 'Q (charge)', 'value': 'q'},
                                    {'label': 'Topology', 'value': 'topology'},
                                    {'label': 'Physics', 'value': 'physics'},
                                ],
                                value='q',
                                style={'color': "#000000", 'width': '50%'}
                            ),
                        ], style={'display': 'flex', 'flexDirection': 'row', 'gap': '10px', 'width': '100%'}),
                        html.H2(),
                        html.Label('Segment/Hit Scale'),
                        dcc.Slider(
                            id='right-window-scale-slider',
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
        ])

    def construct_widgets(self):
        # Callbacks to update the content based on which dropdown is selected
        @self.app.callback(
            Output("dynamic-left-content", "children"),
            [Input("left-window-dropdown", "value")]
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
            Output("dynamic-right-content", "children"),
            [Input("right-window-dropdown", "value")]
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
             Output('tpc_plot_right', 'figure')],
            [Input('event_dropdown', 'value'),
             Input('left-window-plottype-dropdown', 'value'),
             Input('right-window-plottype-dropdown', 'value')],
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
                        charge_segments = flow_file['mc_truth/calib_final_hit_backtrack/data']['segment_ids'].astype(int)
                        charge_fraction = flow_file['mc_truth/calib_final_hit_backtrack/data']['fraction']
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
                        self.charge = flow_file['charge/calib_final_hits/data'][
                            hits_to_segments
                        ]
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
                        arrakis_event_ids = arrakis_file["charge/calib_final_hits/data"]["event_id"]
                        self.topology = arrakis_file["charge/calib_final_hits/data"]["topology"][
                            (arrakis_event_ids == event)
                        ]
                        self.physics = arrakis_file["charge/calib_final_hits/data"]["physics"][
                            (arrakis_event_ids == event)
                        ]
                    self.left_tpc.update_arrakis_event(
                        self.topology,
                        self.physics
                    )
                    self.right_tpc.update_arrakis_event(
                        self.topology,
                        self.physics
                    )
            self.left_tpc.plot_event()
            self.right_tpc.plot_event()
            return print_output, self.left_tpc.tpc, self.right_tpc.tpc

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

        # @self.app.callback(
        #     Output('tpc_plot_left', 'figure'),
        #     Input('left-window-scale-slider', 'value')
        # )
        # def update_left_size_scaler(size):
        #     self.left_tpc.scale = size
        #     print(size)
        #     if self.event is not None:
        #         self.left_tpc.plot_event()
        #     return self.left_tpc.tpc

        # @self.app.callback(
        #     Output('tpc_plot_right', 'figure'),
        #     Input('right-window-scale-slider', 'value')
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
