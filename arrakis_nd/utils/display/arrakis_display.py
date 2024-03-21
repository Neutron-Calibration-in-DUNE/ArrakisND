import os, dash, yaml
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import numpy as np
import h5py
import glob

from arrakis_nd.utils.display.set_server import SetServer
from arrakis_nd.utils.display.vis_event import VisEvent
# from arrakis_nd.utils.display.marjolein_2x2_display import get_layout,add_callbacks

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
        self.arrakis_files = []
        self.available_events = []
        self.event = {"label":-1, "id":-1}
        self.unique_events = []

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

        with open('arrakis_nd/utils/display/styles.yaml', 'r') as file: 
            styles = yaml.safe_load(file)

        # Define the navbar with a dropdown
        # self.navbar = dbc.NavbarSimple(
        #     children=[
        #         dbc.DropdownMenu(
        #             children=[
        #                 dbc.DropdownMenuItem("Item 1", href="#"),
        #                 dbc.DropdownMenuItem("Item 2", href="#"),
        #                 dbc.DropdownMenuItem(divider=True),
        #                 dbc.DropdownMenuItem("Separated Item", href="#"),
        #             ],
        #             nav=True,
        #             in_navbar=True,
        #             label="Menu",
        #         ),
        #     ],
        #     brand="ARRAKIS Display",
        #     brand_href="#",
        #     color="primary",
        #     dark=True,
        # )

        # Define the sidebar
        self.sidebar = html.Div(
            [
                html.H3("ARRAKIS DISPLAY"),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
                
                html.P("🔍 FLOW/ARRAKIS Folders & Files "),
                html.Hr(style={'border': '3px solid #ffffff', 'height': '0px'}),
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
                
                html.Label('FLOW files'),
                dcc.Dropdown(
                    id='flow_dropdown',
                    style={'color': "#000000"}
                ),
                html.Label('ARRAKIS files'),
                dcc.Dropdown(
                    id='arrakis_dropdown',
                    style={'color': "#000000"}
                ),

                html.H2(),
                html.H2(),
                html.Label('Event (spill)'),
                dcc.Dropdown(
                    id='event_dropdown',
                    options=self.available_events,
                    style={'color': "#000000"}
                ),
                html.H2(),
                html.Div(id='event_output'),

                html.H2(),
                html.Hr(),
                html.P(
                    "CHOOSE YOUR DISPLAY", className="lead"
                ),
                
                dbc.Nav(
                    [
                        dbc.NavLink("Event", href="/tab-1", id="tab-1-link", active="exact",style={'color': 'white'}),
                        dbc.NavLink("3d & tree", href="/tab-2", id="tab-2-link", active="exact",style={'color': 'white'}),
                        dbc.NavLink("Tab 3", href="/tab-3", id="tab-3-link", active="exact",style={'color': 'white'}),
                    ],
                    vertical=True,
                    pills=True,
                ),
            ],
            style=styles['SIDEBAR_STYLE'],
        )

        # Define content for the tabs
        self.content = html.Div(id="page-content", style=styles["CONTENT_STYLE"])

        # Define the layout
        self.app.layout = html.Div(style={'overflow': 'scroll'}, children=[
            dcc.Location(id="url"),
            # self.navbar,
            self.sidebar,
            self.content,
        ])

    def construct_widgets(self):
        # Callbacks to update the content based on which tab is selected
        @self.app.callback(
            Output("page-content", "children"),
            [Input("url", "pathname")]
        )
        def render_tab_content(pathname):
            if pathname == "/":
                return html.P("This is the content of the home page!")
            if pathname == "/tab-1":
                # evt_app = DashProxy(__name__, title="2x2 event display")
                # evt_app.layout = get_layout()
                # add_callbacks(evt_app)
                # return get_layout()
                vis_event = VisEvent(self.flow_folder,self.flow_file,self.event)
                return vis_event.get_layout()
                # return create_vis_event()
            elif pathname == "/tab-2":
                return html.P("This is the content of page 2. Yay!")
            elif pathname == "/tab-3":
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
        #TODO: ADD displays stored in vis_*.py files to the tabs

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
            Output('event_dropdown','options'),
            [Input('flow_dropdown', 'value'), Input('arrakis_dropdown', 'value')],
        )
        def update_available_events(flow_file, arrakis_file):
            self.available_events = []
            if flow_file is not None:
                try:
                    self.flow_file = flow_file
                    flow_file = h5py.File(self.flow_folder + flow_file, "r")
                    # self.flow_file = h5flow.data.H5FlowDataManager(self.flow_folder + flow_file, "r")
                    trajectories = flow_file['mc_truth/trajectories/data']
                    events = trajectories['event_id']
                    self.unique_events = np.unique(events)
                    self.available_events = [
                        {'label': event, 'value': event}
                        for event in self.unique_events
                    ]
                except Exception:
                    pass
            return self.available_events

        @self.app.callback(
            Output('event_output', 'children'),
            [Input('event_dropdown', 'value')],
        )
        def load_event(event_file):
            print_output = ''
            if event_file is not None:
                print(f"Event: {event_file}")
                print_output = f"Event: {event_file} Loaded!"
                self.event = {"label":event_file, "id":int(np.where(self.unique_events == event_file)[0])}
            return print_output

    def run_app(self):
        self.app.run_server(
            jupyter_mode=self.jupyter_mode,
            jupyter_server_url=self.server_url,
            host="localhost",
            port=self.port,
        )
        if self.jupyter_mode == "inline":
            self.adjust_iframe_height(height=1500)