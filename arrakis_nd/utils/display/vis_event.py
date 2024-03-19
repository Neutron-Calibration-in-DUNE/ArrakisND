import plotly.graph_objects as go
from dash import dcc,html
import numpy as np
from arrakis_nd.utils.display.utils import *
import h5flow
class VisEvent:
    '''
    This class is used to visualize the events in three different plots:

    '''
    def __init__(self, file, event):
        print("Creating a VisEvent object")
        self.file = file
        self.event = event
        self.data = h5flow.data.H5FlowDataManager(self.file, "r")
        self.num_events = self.data["charge/events/data"].shape[0]
        print(self.data)
        print("Data has been set")

    def create_3d_figure(self):
        fig = go.Figure()
        print("here we go")
        # Select the hits for the current event
        prompthits_ev = self.data["charge/events", "charge/calib_prompt_hits", self.event["id"]]
        finalhits_ev = self.data["charge/events", "charge/calib_final_hits", self.event["id"]]
        # select the segments (truth) for the current event
        try:
            prompthits_segs = self.data[
                "charge/events",
                "charge/calib_prompt_hits",
                "charge/packets",
                "mc_truth/segments",  # called segments in minirun4
                self.event["id"],
            ]
            sim_version = "minirun4"
            print("Found truth info in minirun4 format")
        except:
            print("No truth info in minirun4 format found")
            try:
                prompthits_segs = self.data[
                    "charge/events",
                    "charge/calib_prompt_hits",
                    "charge/packets",
                    "mc_truth/tracks",  # called tracks in minirun3
                    self.event["id"],
                ]
                sim_version = "minirun3"
                print("Found truth info in minirun3 format")
            except:
                print("No truth info in minirun3 format found")
                prompthits_segs = None

        # Plot the prompt hits
        print("Plotting prompt hits")
        print(prompthits_ev)
        print(prompthits_ev["x"])
        print(prompthits_ev.data["x"])
        prompthits_traces = go.Scatter3d(
            x=prompthits_ev.data["x"].flatten(),
            y=(prompthits_ev.data["y"].flatten()),
            z=(prompthits_ev.data["z"].flatten()),
            marker_color=prompthits_ev.data["E"].flatten()
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
            customdata=prompthits_ev.data["E"].flatten() * 1000,
            hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
        )
        print("Adding prompt hits to figure")
        fig.add_traces(prompthits_traces)

        # Plot the final hits
        print("Plotting final hits")
        finalhits_traces = go.Scatter3d(
            x=finalhits_ev.data["x"].flatten(),
            y=(finalhits_ev.data["y"].flatten()),
            z=(finalhits_ev.data["z"].flatten()),
            marker_color=finalhits_ev.data["E"].flatten() * 1000,
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
            customdata=finalhits_ev.data["E"].flatten() * 1000,
            hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
        )
        print("Adding final hits to figure")
        fig.add_traces(finalhits_traces)

        print("Plotting segments")
        if prompthits_segs is not None:
            segs_traces = plot_segs(
                prompthits_segs[0, :, 0, 0],
                sim_version=sim_version,
                mode="lines",
                name="edep segments",
                visible="legendonly",
                line_color="red",
                showlegend=True,
            )
            print("Adding segments to figure")
            fig.add_traces(segs_traces)

        # Draw the TPC
        print("Drawing TPC")
        tpc_center, anodes, cathodes = draw_tpc(sim_version)
        light_detectors = draw_light_detectors(self.file, self.file["id"])

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

### Old code (works)
    # def get_layout(self):
    #     print("Creating the layout for the VisEvent object")
    #     # Create a 3D scatter plot with self.data
    #     fig = go.Figure(data=[go.Scatter3d(
    #         x=np.array(self.data['x']),
    #         y=np.array(self.data['y']),
    #         z=np.array(self.data['z']),
    #         mode='markers',
    #         marker=dict(
    #             size=12,
    #             # color=self.data['color'],  # set color to an array/list of desired values
    #             colorscale='Viridis',  # choose a colorscale
    #             opacity=0.8
    #         )
    #     )])

    #     # Set the layout of the plot
    #     fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))

    #     # Return a Div containing the plot
    #     return html.Div([
    #         dcc.Graph(
    #             id='3d-scatter',
    #             figure=fig
    #         )
    #     ])
