import os, dash, yaml
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

from arrakis_nd.utils.display.set_server import SetServer
from arrakis_nd.utils.display.vis_event import create_3d_scatter

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

        

        self.available_events = []
        self.event = -1

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
        self.sidebar = html.Div([html.H3("ARRAKIS DISPLAY"),
                html.Hr(),
                html.P("🔍 Folders & Files "),
                html.Hr(),
                dbc.Label("flow_folder"),
                dbc.Input(placeholder="Enter your flow_folder",
                            type="text", id='flow_folder_input'),
                dbc.Label("arrakis_folder"),
                dbc.Input(placeholder="Enter your arrakis_folder",
                            type="text", id='arrakis_folder_input'),

                html.H2(),
                dbc.Button("Load files", color="primary", 
                                className="me-1", id='update_files'),
                html.H2(),
                html.Label('flow_file'),
                dcc.Dropdown(id='flow_dropdown'),
                html.Label('arrakis_file'),
                dcc.Dropdown(id='arrakis_dropdown'),

                html.H2(),
                dbc.Label("Choose one"),
                dbc.RadioItems(
                    options=
                    [
                        {"label": "Flow", "value": "flow_file"},
                        {"label": "Arrakis", "value": "arrakis_file"},
                    ],
                    value=None,
                    id="radioitems_file",
                    inline=True,
                ),
                dcc.Dropdown(id='event_dropdown',options=self.available_events),
                html.H2(),
                dbc.Button("Load event", color="primary", 
                        className="me-1", id='update_event'),
                html.Div(id='event_output'),

                html.H2(),
                html.Hr(),
                html.P(
                    "CHOOSE YOUR DISPLAY", className="lead"
                ),
                
                dbc.Nav(
                    [
                        dbc.NavLink("Event", href="/tab-1", id="tab-1-link", active="exact"),
                        dbc.NavLink("3d & tree", href="/tab-2", id="tab-2-link", active="exact"),
                        dbc.NavLink("Tab 3", href="/tab-3", id="tab-3-link", active="exact"),
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

    def load_event(self, event):
        pass

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
                return create_3d_scatter()
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
            [Output('flow_dropdown', 'options'), Output('arrakis_dropdown', 'options')],
            [Input('update_files', 'n_clicks')],
            [dash.dependencies.State('flow_folder_input', 'value'),
            dash.dependencies.State('arrakis_folder_input', 'value')]
        )
        def update_dropdown(n_clicks, flow_folder, arrakis_folder):
            if n_clicks is not None:
                flow_options = []
                arrakis_options = []
                if flow_folder and os.path.isdir(flow_folder):
                    flow_files = os.listdir(flow_folder)
                    # print(f"Flow Files: {flow_files}")  # Print the list of flow files
                    flow_options = [{'label': file, 'value': file} for file in flow_files]
                if arrakis_folder and os.path.isdir(arrakis_folder):
                    arrakis_files = os.listdir(arrakis_folder)
                    # print(f"Arrakis Files: {arrakis_files}")  # Print the list of arrakis files
                    arrakis_options = [{'label': file, 'value': file} for file in arrakis_files]
                return flow_options, arrakis_options
            return [], []
        
        @self.app.callback(
            Output('event_dropdown','options'),
            [Input('radioitems_file', 'value')],
            [dash.dependencies.State('flow_dropdown', 'value'),
            dash.dependencies.State('arrakis_dropdown', 'value')]
        )
        def update_dropdown(radio_value, flow_file, arrakis_file):
            self.available_events = []
            print(f"Radio Value: {radio_value}")  # Print the selected radio item
            if radio_value == 'flow_file':
                print(f"Flow File: {flow_file}")  # Print the selected flow file
                # Load events from the selected flow file
                self.available_events = [{'label': 'Event 1', 'value': 1}, {'label': 'Event 2', 'value': 2}]
            elif radio_value == 'arrakis_file':
                print(f"Arrakis File: {arrakis_file}")  # Print the selected arrakis file
                # Load events from the selected arrakis file
                self.available_events = [{'label': 'Event 1', 'value': 1}, {'label': 'Event 2', 'value': 2}]
            return self.available_events
        #TODO: UPDATE WITH AN EVENT LOADER FUNCTION TO GET THE EVENTS FROM THE FILES

        @self.app.callback(
            Output('event_output', 'children'),
            [Input('update_event', 'n_clicks')],
            dash.dependencies.State('event_dropdown', 'value'),
        )
        def update_dropdown(n_clicks, event_file):
            print_output = ''
            if n_clicks is not None:
                print(f"Event: {event_file}")
                print_output = f"Event: {event_file} Loaded!"
                # self.event = event_file
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