import os, socket, dash
from dash import Dash, dcc, html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np


class ArrakisDisplay:
    """
    A class to create a Dash app for displaying the results of the ArrakisND analysis.
    """
    def __init__(
        self,
        service_prefix: str = os.getenv("JUPYTERHUB_SERVICE_PREFIX", "/"),
        server_url: str = "https://jupyter.nersc.gov",
        port: int = 8050,
    ):
        """
        Initialize the ArrakisDisplay class.

        Args:
            service_prefix (str, optional): _description_. Defaults to os.getenv("JUPYTERHUB_SERVICE_PREFIX", "/").
            server_url (_type_, optional): _description_. Defaults to "https://jupyter.nersc.gov".
            port (int, optional): _description_. Defaults to 8050.
        """
        self.service_prefix = service_prefix
        self.server_url = server_url
        self.port = self.find_free_port(start_port=port)

        self.available_events = []
        self.event = -1

        self.construct_app()
        self.construct_widgets()
        self.run_app()

    def find_free_port(
        self,
        start_port=8050
    ):
        """
        Find a free port to run the Dash app on.

        Args:
            start_port (int, optional): _description_. Defaults to 8050.

        Raises:
            OSError: _description_

        Returns:
            _type_: _description_
        """
        port = start_port
        while True:
            try:
                sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                sock.bind(("localhost", port))
                sock.close()
                return port
            except OSError:
                port += 1
                if port - start_port > 1000:  # Just to avoid an infinite loop
                    raise OSError("No free ports available within 1000 ports of start_port")

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
        self.app = Dash(
            __name__,
            requests_pathname_prefix=f"{self.service_prefix}proxy/{self.port}/",
            external_stylesheets=[dbc.themes.SPACELAB]
        )
        # Define the navbar with a dropdown
        self.navbar = dbc.NavbarSimple(
            children=[
                dbc.DropdownMenu(
                    children=[
                        dbc.DropdownMenuItem("Item 1", href="#"),
                        dbc.DropdownMenuItem("Item 2", href="#"),
                        dbc.DropdownMenuItem(divider=True),
                        dbc.DropdownMenuItem("Separated Item", href="#"),
                    ],
                    nav=True,
                    in_navbar=True,
                    label="Menu",
                ),
            ],
            brand="ARRAKIS Display",
            brand_href="#",
            color="primary",
            dark=True,
        )

        # Define the sidebar
        self.sidebar = html.Div(
            [
                html.H2("Sidebar", className="display-4"),
                html.Hr(),
                html.P(
                    "Introduce the folders with the data you want to visualize", className="lead"
                ),
                html.Div(
                    [
                        dbc.Label("flow_folder"),
                        dbc.Input(placeholder="Enter your flow_folder", type="text", id='flow_folder_input'),
                        # dbc.FormText("Type something in the box above"),
                    ]
                ),

                html.Div(
                    [
                        dbc.Label("arrakis_folder"),
                        dbc.Input(placeholder="Enter your arrakis_folder", type="text", id='arrakis_folder_input'),
                        # dbc.FormText("Type something in the box above"),
                    ]
                ),

                html.H3(),
                html.Div(
                    [   
                        dbc.Button("Load files", color="primary", className="me-1", id='update_files'),
                    ]
                ),
                html.H3(),

                html.Label('flow_file'),
                dcc.Dropdown(id='flow_dropdown'),
                html.Label('arrakis_file'),
                dcc.Dropdown(id='arrakis_dropdown'),

                html.H3(),
                html.Div(
                    [
                        dbc.Label("Choose one"),
                        dbc.RadioItems(
                            options=[
                                {"label": "Flow", "value": "flow_file"},
                                {"label": "Arrakis", "value": "arrakis_file"},
                            ],
                            value=None,
                            id="radioitems_file",
                            inline=True,
                        ),
                    ]
                ),

                html.Div(
                    [   
                        dcc.Dropdown(id='event_dropdown',options=self.available_events),
                        dbc.Button("Load event", color="primary", className="me-1", id='update_event'),
                        html.Div(id='event_output')
                    ]
                ),

                html.H3(),
                html.Hr(),
                html.P(
                    "Choose the display you want to show", className="lead"
                ),
                dbc.Nav(
                    [
                        dbc.NavLink("Event", href="/tab-1", id="tab-1-link"),
                        dbc.NavLink("3d & tree", href="/tab-2", id="tab-2-link"),
                        dbc.NavLink("Tab 3", href="/tab-3", id="tab-3-link"),
                    ],
                    vertical=True,
                    pills=True,
                ),
            ],
            style={"width": "20%", "height": "100vh", "position": "fixed", "padding": "2rem 1rem"},
        )

        # Define content for the tabs
        self.tab_content = html.Div(
            id="tab-content",
            children=[
                dcc.Tabs(id="tabs", children=[
                    dcc.Tab(label="Tab 1", children=[html.Div("EVENT DISPLAY")], id="tab-1"),
                    dcc.Tab(label="Tab 2", children=[html.Div("3D & TREE CLUSTERING DISPLAY")], id="tab-2"),
                    dcc.Tab(label="Tab 3", children=[html.Div("Content for Tab 3")], id="tab-3"),
                ], style={"margin-left": "20%"}),
            ]
        )

        # Define the layout
        self.app.layout = html.Div([
            self.navbar,
            self.sidebar,
            self.tab_content,
        ])

    def load_event(self, event):
        pass

    def construct_widgets(self):
        # Callbacks to update the content based on which tab is selected
        @self.app.callback(
            Output("tab-content", "children"),
            [
                Input("tab-1-link", "n_clicks"),
                Input("tab-2-link", "n_clicks"),
                Input("tab-3-link", "n_clicks")
            ]
        )
        def render_tab_content(tab1, tab2, tab3):
            ctx = dash.callback_context

            if not ctx.triggered:
                return "EVENT DISPLAY"
            else:
                button_id = ctx.triggered[0]["prop_id"].split(".")[0]

                if button_id == "tab-1-link":
                    return "EVENT DISPLAY"
                elif button_id == "tab-2-link":
                    return "3D & TREE CLUSTERING DISPLAY"
                elif button_id == "tab-3-link":
                    return "This is Tab 3's content"

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
            jupyter_mode="inline",
            jupyter_server_url=self.server_url,
            host="localhost",
            port=self.port
        )
        self.adjust_iframe_height(height=1500)
