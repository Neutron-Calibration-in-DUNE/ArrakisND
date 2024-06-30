import plotly.graph_objects as go
from dash import dcc,html
from arrakis_nd.utils.display.utils import *
import h5py
class VisEvent:
    '''
    This class is used to visualize the events in three different plots:

    '''
    def __init__(self,folder,file,event):
        print("Creating a VisEvent object")

        self.folder = folder
        self.file = file
        self.event = event
        self.data = h5py.File(self.folder + file, "r")
        self.sim_version = ''
        

    def create_3d_figure(self):
        fig = go.Figure()
        print("here we go")
        # Select the hits for the current event
        prompthits_ev = self.data["charge/calib_prompt_hits/data"]         
        finalhits_ev = self.data["charge/calib_final_hits/data"]
        # select the segments (truth) for the current event
        try:
            prompthits_segs = self.data["mc_truth/segments"]
            self.sim_version = "minirun4"
            print("Found truth info in minirun4 format")
        except:
            print("No truth info in minirun4 format found")
            try:
                prompthits_segs = self.data["mc_truth/tracks"]
                self.sim_version = "minirun3"
                print("Found truth info in minirun3 format")
            except:
                print("No truth info in minirun3 format found")
                prompthits_segs = None
        # TODO: understand/fix this
                
        # Plot the prompt hits
        print("Plotting prompt hits")
        prompthits_traces = go.Scatter3d(
            x=prompthits_ev["x"].flatten(),
            y=(prompthits_ev["y"].flatten()),
            z=(prompthits_ev["z"].flatten()),
            marker_color=prompthits_ev["E"].flatten()
            * 1000,  # convert to MeV from GeV for minirun4, not sure for minirun3
            marker={
                "size": 1.75,
                "opacity": 0.7,
                "colorscale": "cividis",
                "colorbar": {
                    "title": "Hit energy [MeV]",
                    "titlefont": {"size": 12},
                    "tickfont": {"size": 10},
                    "thickness": 15,
                    "len": 0.5,
                    "xanchor": "left",
                    "x": 0,
                },
            },
            name="prompt hits",
            mode="markers",
            showlegend=True,
            opacity=0.7,
            customdata=prompthits_ev["E"].flatten() * 1000,
            hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
        )
        print("Adding prompt hits to figure")
        fig.add_traces(prompthits_traces)

        # Plot the final hits
        print("Plotting final hits")
        finalhits_traces = go.Scatter3d(
            x=finalhits_ev["x"].flatten(),
            y=(finalhits_ev["y"].flatten()),
            z=(finalhits_ev["z"].flatten()),
            marker_color=finalhits_ev["E"].flatten() * 1000,
            marker={
                "size": 1.75,
                "opacity": 0.7,
                "colorscale": "Plasma",
                "colorbar": {
                    "title": "Hit energy [MeV]",
                    "titlefont": {"size": 12},
                    "tickfont": {"size": 10},
                    "thickness": 15,
                    "len": 0.5,
                    "xanchor": "left",
                    "x": 0,
                },
            },
            name="final hits",
            mode="markers",
            visible="legendonly",
            showlegend=True,
            opacity=0.7,
            customdata=finalhits_ev["E"].flatten() * 1000,
            hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
        )
        print("Adding final hits to figure")
        fig.add_traces(finalhits_traces)

        # print("Plotting segments")
        # if prompthits_segs is not None:
        #     segs_traces = plot_segs(
        #         # prompthits_segs[0, :, 0, 0],
        #         sim_version=self.sim_version,
        #         mode="lines",
        #         name="edep segments",
        #         visible="legendonly",
        #         line_color="red",
        #         showlegend=True,
        #     )
        #     print("Adding segments to figure")
        #     fig.add_traces(segs_traces)
        # TODO: fix this

        # Draw the TPC
        print("Drawing TPC")
        tpc_center, anodes, cathodes = draw_tpc(self.sim_version)
        light_detectors = draw_light_detectors(self.data, 0)

        print("Adding TPC to figure")
        fig.add_traces(tpc_center)
        fig.add_traces(anodes)
        fig.add_traces(cathodes)
        print("Adding light detectors to figure")
        fig.add_traces(light_detectors)

        return fig

    def get_layout(self):
        figure3d = self.create_3d_figure()
        return html.Div(
            [
                # Hidden divs to store data
                dcc.Location(id="url"),
                dcc.Store(id="filename", storage_type="local", data=None),
                dcc.Store(id='data-length', data=0),
                # Header
                html.H1(children="2x2 event display", style={"textAlign": "center"}),
                # Graphs
                html.Div([
                    # Existing 3D graph
                    html.Div(dcc.Graph(id='3d-graph', style={'height': '70vh', 'width': '50vw'}, figure=figure3d)),
                    # New Light waveform graph
                    html.Div(dcc.Graph(id="light-waveform", style={'height': '50vh', 'width': '35vw'})),
                ], style={'display': 'flex'}),
                # New Another Graph (replace with your actual component)
                html.Div(dcc.Graph(id="another-graph", style={'height': '30vh', 'width': '35vw', 'float': 'right'})),
            ]
        )