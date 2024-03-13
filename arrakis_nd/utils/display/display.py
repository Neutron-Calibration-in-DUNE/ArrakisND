import os
import socket
import dash
from dash import Dash
from dash import dcc
from dash import html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc
import plotly.graph_objs as go
import numpy as np


class ArrakisDisplay:
    """_summary_
    """
    def __init__(
        self,
        service_prefix: str = os.getenv("JUPYTERHUB_SERVICE_PREFIX", "/"),
        server_url: str = "https://jupyter.nersc.gov",
        port: int = 8050,
    ):
        """_summary_

        Args:
            service_prefix (str, optional): _description_. Defaults to os.getenv("JUPYTERHUB_SERVICE_PREFIX", "/").
            server_url (_type_, optional): _description_. Defaults to "https://jupyter.nersc.gov".
            port (int, optional): _description_. Defaults to 8050.
        """
        self.service_prefix = service_prefix
        self.server_url = server_url
        self.port = self.find_free_port(start_port=port)

        self.construct_app()
        self.construct_widgets()
        self.run_app()

    def find_free_port(
        self,
        start_port=8050
    ):
        """_summary_

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
        """_summary_
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
                    "A simple sidebar layout with navigation links", className="lead"
                ),
                dbc.Nav(
                    [
                        dbc.NavLink("Tab 1", href="/tab-1", id="tab-1-link"),
                        dbc.NavLink("Tab 2", href="/tab-2", id="tab-2-link"),
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
                    dcc.Tab(label="Tab 1", children=[html.Div("Content for Tab 1")], id="tab-1"),
                    dcc.Tab(label="Tab 2", children=[html.Div("Content for Tab 2")], id="tab-2"),
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
                return "This is Tab 1's content"
            else:
                button_id = ctx.triggered[0]["prop_id"].split(".")[0]

                if button_id == "tab-1-link":
                    return "This is Tab 1's content"
                elif button_id == "tab-2-link":
                    return "This is Tab 2's content"
                elif button_id == "tab-3-link":
                    return "This is Tab 3's content"

    def run_app(self):
        self.app.run_server(
            jupyter_mode="inline",
            jupyter_server_url=self.server_url,
            host="localhost",
            port=self.port
        )
        self.adjust_iframe_height(height=1500)
