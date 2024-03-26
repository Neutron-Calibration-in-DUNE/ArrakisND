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
        self.interactions = None
        self.segments = None
        self.stack = None
        self.trajectories = None
        self.charge = None
        self.set_geometry_info(geometry_info)
        self.generate_layout()
    
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
        self.tpc = self.construct_tpcs()
        self.waveforms = self.construct_waveforms()
        self.larpix = self.construct_larpix()
    
    def construct_tpcs(self):
        tpc = go.Figure()
        tpc.update_layout(
            scene=dict(
                xaxis_title="x [cm]",
                yaxis_title="z [cm]",
                zaxis_title="y [cm]"
            ),
            title="2x2 TPCs"
        )
        self.tpc_center, self.anodes, self.cathodes = self.draw_tpc()
        tpc.add_traces(self.tpc_center)
        tpc.add_traces(self.anodes)
        tpc.add_traces(self.cathodes)
        return tpc
    
    def construct_light_detectors(self):
        pass
    
    def construct_waveforms(self):
        waveforms = go.Figure()
        return waveforms
        
    def construct_larpix(self):
        larpix = go.Figure()
        return larpix
    
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
                traces.append(go.Surface(x=x, y=z, z=y, **kwargs))

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

                traces.append(go.Surface(x=x, y=z, z=y, **kwargs))

        return traces
        
    def update_event(
        self,
        interactions,
        segments,
        stack,
        trajectories,
        charge
    ):
        self.interactions = interactions
        self.segments = segments
        self.stack = stack
        self.trajectories = trajectories
        self.charge = charge
    
    def plot_event(self):
        if self.charge is not None:
            name_to_remove = "calib_final_hits"
            # Create a new list of traces that excludes the ones with the specified name
            new_traces = [trace for trace in self.tpc.data if trace.name != name_to_remove]
            # Update the figure with the new list of traces
            self.tpc.data = new_traces
            
            charge_hits_traces = go.Scatter3d(
                x=self.charge['x'],
                y=self.charge['z'],
                z=self.charge['y'],
                marker_color=self.charge['Q'],
                marker={
                    "size": 1.75,
                    "opacity": 0.7,
                    "colorscale": "cividis",
                    "colorbar": {
                        "title": "Hit charge [e-]",
                        "titlefont": {"size": 12},
                        "tickfont": {"size": 10},
                        "thickness": 15,
                        "len": 0.5,
                        "xanchor": "left",
                        "x": 0,
                    },
                },
                name="calib_final_hits",
                mode="markers",
                showlegend=True,
                opacity=0.7,
                hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
            )
            
            self.tpc.add_traces(charge_hits_traces)

    def generate_layout(self):
        self.construct_detector()
        self.plot_event()
        self.layout = html.Div(
            [
                dcc.Location(id="url"),
                dcc.Store(id="filename", storage_type="local", data=None),
                dcc.Store(id='data-length', data=0),
                html.H1(children="2x2 event display", style={"textAlign": "center"}),
                html.Div([
                    html.Div(dcc.Graph(
                        id='charge_light_tpc_plot', 
                        style={'height': '70vh', 'width': '50vw'}, 
                        figure=self.tpc)
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