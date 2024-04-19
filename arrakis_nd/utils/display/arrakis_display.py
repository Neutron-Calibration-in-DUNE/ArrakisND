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
from arrakis_nd.utils.display.charge_light_display import ChargeLightDisplay
from arrakis_nd.utils.display.vis_event import VisEvent


class ArrakisDisplay:
    """
    A class to create a Dash app for displaying the results of the ArrakisND analysis.
    """
    def __init__(
        self,
    ):
        """
        Initialize the ArrakisDisplay class.

        Args:
            service_prefix (str, optional): _description_. Defaults to os.getenv("JUPYTERHUB_SERVICE_PREFIX", "/").
            server_url (_type_, optional): _description_. Defaults to "https://jupyter.nersc.gov".
            port (int, optional): _description_. Defaults to 8050.
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

        self.standard_flow_folders = [
            {'label': 'MiniRun4', 'value':
                '/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun4_1E19_RHC/MiniRun4_1E19_RHC.flow/FLOW'},
            {'label': 'MiniRun4.5', 'value':
                '/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun4.5_1E19_RHC/inputs_beta3/MiniRun4.5_1E19_RHC.flow/FLOW/0000000'},
            {'label': 'MiniRun5', 'value':
                '/global/cfs/cdirs/dune/www/data/2x2/simulation/productions/MiniRun5_1E19_RHC/MiniRun5_1E19_RHC.flow.beta1/FLOW/0000000'},
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

        self.charge_light_display = ChargeLightDisplay()

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

                # Event display selector
                html.H2(),
                html.Hr(),
                html.P(
                    "Choose a Display",
                    className="lead"
                ),
                dcc.Dropdown(
                    id='display_dropdown',
                    options=[
                        'Home',
                        'Charge and Light Detectors',
                        'Arrakis',
                        'Clustering'
                    ],
                    value='Home',
                    style={'color': "#000000"}
                ),

                # Highlighting options
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                html.P("FLOW charge types to plot"),
                dcc.Checklist(
                    ['Segments', 'Calib Prompt Hits', 'Calib Final Hits'],
                    ['Calib Prompt Hits'], inline=True
                ),
                html.Hr(style={'border': '1px solid #ffffff', 'height': '0px'}),
                html.Div([
                    html.Div([
                        html.P("Active TPCs"),
                        dcc.Checklist(
                            ['0', '1', '2', '3', '4', '5', '6', '7'],
                            ['0', '1', '2', '3', '4', '5', '6', '7']
                        )], style={'width': '45vw'}),
                    html.Div([
                        html.P("Detectors"),
                        dcc.Checklist(
                            ['0', '1', '2', '3', '4', '5', '6', '7'],
                            ['0', '1', '2', '3', '4', '5', '6', '7']
                        )], style={'width': '45vw'}),
                ], style={'display': 'flex'}),
            ],
            style=styles['SIDEBAR_STYLE'],
        )

        # Define content for the tabs
        self.content = html.Div(id="page-content", style=styles["CONTENT_STYLE"])

        # Define the layout
        self.app.layout = html.Div(style={'overflow': 'scroll'}, children=[
            dcc.Location(id="url"),
            self.navbar,
            self.sidebar,
            self.content,
        ])

    def construct_widgets(self):
        # Callbacks to update the content based on which tab is selected
        @self.app.callback(
            Output("page-content", "children"),
            [Input("display_dropdown", "value")]
        )
        def render_tab_content(pathname):
            if pathname == "Home":
                return html.P("This is the content of the home page!")
            if pathname == "Charge and Light Detectors":
                print("Charge and Light Detectors")
                # print(self.charge_light_display.layout)
                return self.charge_light_display.layout
            elif pathname == "Arrakis":
                return html.P("This is the content of page 2. Yay!")
            elif pathname == "Clustering":
                return html.P("Oh cool, this is page 3!")
            # If the user tries to reach a different page, return a 404 message
            return html.Div(
                [
                    html.H1("404: Not found", className="text-danger"),
                    html.Hr(),
                    html.P(f"The pathname {pathname} was not recognised..."),
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
                        self.charge_light_display.set_geometry_info(self.geometry_info)
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
            return self.available_events[new_index]['value']

        @self.app.callback(
            [Output('event_output', 'children'),
             Output('charge_light_tpc_plot', 'figure')],
            [Input('event_dropdown', 'value')],
        )
        def load_event(event):
            print_output = ''
            if event is not None:
                if self.flow_file:
                    with h5py.File(self.flow_folder + self.flow_file, "r") as flow_file:
                        interactions_events = flow_file['mc_truth/interactions/data']['event_id']
                        #print(np.where(interactions_events == event)[0])
                        segments_events = flow_file['mc_truth/segments/data']['event_id']
                        stack_events = flow_file['mc_truth/stack/data']['event_id']
                        trajectories_events = flow_file['mc_truth/trajectories/data']['event_id']
                        charge_segments = flow_file['mc_truth/calib_final_hit_backtrack/data']['segment_id'].astype(int)
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
                        charge_events = flow_file["charge/events/data"]["id"]

                        self.charge_events = flow_file["charge/events/data"][np.where(interactions_events == event)[0]]
                        
                        """Likewise for light data, we must backtrack through segments"""
                        match_light = flow_file['/light/events/data'][:][ flow_file['/charge/events/ref/light/events/ref'][np.where(interactions_events == event)[0],1] ]["id"]
                        self.light = flow_file["light/events/data"][:]  # we have to try them all, events may not be time ordered

                        waveforms_all_detectors = flow_file["light/wvfm/data"]["samples"][match_light]
                        print(match_light)
                        # we have now the waveforms for all detectors matched in time to the event
                    print_output = f"Event: {event} Loaded!"
                    self.charge_light_display.update_event(
                        self.interactions,
                        self.segments,
                        self.stack,
                        self.trajectories,
                        self.charge,
                    )
                    self.charge_light_display.plot_event()
                    self.charge_light_display.construct_light_detectors(waveforms_all_detectors)
                    self.charge_light_display.waveforms = self.charge_light_display.construct_waveforms()
            return print_output, self.charge_light_display.tpc
        
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
                            match_light = flow_file['/light/events/data'][:][ flow_file['/charge/events/ref/light/events/ref'][np.where(interactions_events == event)[0],1] ]["id"]
                            
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

    def run_app(self):
        self.app.run_server(
            jupyter_mode=self.jupyter_mode,
            jupyter_server_url=self.server_url,
            host="localhost",
            port=self.port,
        )
        if self.jupyter_mode == "inline":
            self.adjust_iframe_height(height=1500)
