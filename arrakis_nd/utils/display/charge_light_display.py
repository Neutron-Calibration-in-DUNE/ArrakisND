import plotly.graph_objects as go
from dash import dcc,html
import numpy as np
import plotly
import plotly.graph_objects as go
import plotly.graph_objects as go

# from arrakis_nd.utils.display.utils import *
import h5py
class ChargeLightDisplay:
    '''
    This class is used to visualize the events in three different plots:

    '''
    def __init__(
        self,
        geometry_info: dict={}
    ):
        self.set_geometry_info(geometry_info)
    
    def set_geometry_info(
        self,
        geometry_info,
    ):
        if geometry_info == {}:
            pass
        else:
            try:
                pass
            except Exception:
                pass
    
    def construct_detector(self):
        print("construct_detector")
        tpc = self.construct_tpcs()
        waveforms = self.construct_waveforms()
        larpix = self.construct_larpix()
        return tpc, waveforms, larpix
    
    def construct_tpcs(self):
        tpc = go.Figure()
        tpc_center, anodes, cathodes = self.draw_tpc()
        tpc.add_traces(tpc_center)
        tpc.add_traces(anodes)
        tpc.add_traces(cathodes)
        return tpc
    
    def construct_light_detectors(self):
        pass
    
    def construct_waveforms(self):
        waveforms = go.Figure()
        return waveforms
        
    def construct_larpix(self):
        larpix = go.Figure()
        return larpix
        
    def plot_event(
        self,
    ):
        pass
    
    def draw_tpc(
        self,
        sim_version="minirun5"
    ):
        anode_xs = np.array([-63.931, -3.069, 3.069, 63.931])
        anode_ys = np.array([-19.8543, 103.8543])  # two ys
        anode_zs = np.array([-64.3163, -2.6837, 2.6837, 64.3163])  # four zs
        if sim_version == "minirun4":  # hit coordinates are in cm
            detector_center = (0, -268, 1300)
            anode_ys = anode_ys - (268 + 42)
            anode_zs = anode_zs + 1300
        if sim_version == "minirun3":  # hit coordinates are in mm
            detector_center = (0, 42 * 10, 0)
            anode_xs = anode_xs * 10
            anode_ys = anode_ys * 10
            anode_zs = anode_zs * 10
        if sim_version == "minirun5":  # hit coordinates are in cm
            detector_center = (0, 0, 0)
            anode_ys = anode_ys - 42

        center = go.Scatter3d(
            x=[detector_center[0]],
            y=[detector_center[1]],
            z=[detector_center[2]],
            marker=dict(size=3, color="green", opacity=0.5),
            mode="markers",
            name="tpc center",
        )
        anodes = self.draw_anode_planes(
            anode_xs, anode_ys, anode_zs, colorscale="ice", showscale=False, opacity=0.1
        )
        cathodes = self.draw_cathode_planes(
            anode_xs, anode_ys, anode_zs, colorscale="burg", showscale=False, opacity=0.1
        )
        return center, anodes, cathodes

    def draw_cathode_planes(
        self,
        x_boundaries, 
        y_boundaries, 
        z_boundaries, 
        **kwargs
    ):
        traces = []
        for i_z in range(int(len(z_boundaries) / 2)):
            for i_x in range(int(len(x_boundaries) / 2)):
                z, y = np.meshgrid(
                    np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2),
                    np.linspace(y_boundaries.min(), y_boundaries.max(), 2),
                )
                x = (
                    (x_boundaries[i_x * 2] + x_boundaries[i_x * 2 + 1])
                    * 0.5
                    * np.ones(z.shape)
                )
                traces.append(go.Surface(x=x, y=y, z=z, **kwargs))

        return traces

    def draw_anode_planes(
        self,
        x_boundaries, 
        y_boundaries, 
        z_boundaries, 
        **kwargs
    ):
        traces = []
        for i_z in range(int(len(z_boundaries) / 2)):
            for i_x in range(int(len(x_boundaries))):
                z, y = np.meshgrid(
                    np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2),
                    np.linspace(y_boundaries.min(), y_boundaries.max(), 2),
                )
                x = x_boundaries[i_x] * np.ones(z.shape)

                traces.append(go.Surface(x=x, y=y, z=z, **kwargs))

        return traces
        
        
    def create_3d_figure(self):
        fig = go.Figure()
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
        light_detectors = draw_light_detectors(self.data,self.event["id"])

        print("Adding TPC to figure")
        fig.add_traces(tpc_center)
        fig.add_traces(anodes)
        fig.add_traces(cathodes)
        print("Adding light detectors to figure")
        fig.add_traces(light_detectors)

        return fig

    def get_layout(self):
        print("get_layout")
        tpc, waveforms, larpix = self.construct_detector()
        return html.Div(
            [
                dcc.Location(id="url"),
                dcc.Store(id="filename", storage_type="local", data=None),
                dcc.Store(id='data-length', data=0),
                html.H1(children="2x2 event display", style={"textAlign": "center"}),
                html.Div([
                    html.Div(dcc.Graph(
                        id='3d-graph', 
                        style={'height': '70vh', 'width': '50vw'}, 
                        figure=tpc)
                    ),
                    # html.Div(dcc.Graph(
                    #     id="light-waveform", 
                    #     style={'height': '50vh', 'width': '35vw'},
                    #     figure=self.waveforms)
                    # ),
                ], style={'display': 'flex'}),
                html.Div(dcc.Graph(id="another-graph", style={'height': '30vh', 'width': '35vw', 'float': 'right'})),
            ]
        )