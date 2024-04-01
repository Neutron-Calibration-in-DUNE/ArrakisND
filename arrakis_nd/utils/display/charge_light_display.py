import plotly.graph_objects as go
from dash import dcc,html
import numpy as np
import plotly
import plotly.graph_objects as go
import plotly.graph_objects as go
import h5py


class ChargeLightDisplay:
    '''
    This class is used to visualize the events in three different plots:

    '''
    def __init__(
        self,
        geometry_info: dict = {}
    ):
        self.interactions = None
        self.segments = None
        self.stack = None
        self.trajectories = None
        self.charge = None
        self.geometry_info = {}
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
                self.geometry_info = geometry_info
                self.light_traps = self.construct_light_detectors()
                self.tpc.add_traces(self.light_traps)
            except Exception as e:
                print(f"issue with light traps {e}")

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
        """Plot optical detectors"""
        drawn_objects = []
        ys = np.flip(np.array(
            [
                -595.43, -545.68, -490.48, -440.73,
                -385.53, -335.78, -283.65, -236.65,
                -178.70, -131.70,  -73.75,  -26.75,
                25.38,   75.13,  130.33,  180.08,
                235.28,  285.03,  337.15,  384.15,
                442.10,  489.10,  547.05,  594.05,
            ]) / 10
        )
        light_width = ys[1] - ys[0]

        det_bounds = self.geometry_info["det_bounds"]
        COLORSCALE = plotly.colors.make_colorscale(
            plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.YlOrRd)[0]
        )

        for ix in range(0, det_bounds.shape[0]):
            for ilight, light_y in enumerate(ys):
                for iside in range(2):
                    opid = ilight + iside * len(ys) + ix * len(ys) * 2
                    # if opid not in op_indeces:
                    #     continue
                    xx = np.linspace(det_bounds[ix][0][0][0], det_bounds[ix][0][1][0], 2)
                    zz = np.linspace(
                        light_y - light_width / 2 + det_bounds[0][0][1] + 0.25,
                        light_y + light_width / 2 + det_bounds[0][0][1] - 0.25,
                        2,
                    )
                    xx, zz = np.meshgrid(xx, zz)
                    # light_color = [
                    #     [
                    #         0.0,
                    #         get_continuous_color(
                    #             COLORSCALE, intermed=max(0, n_photons[opid]) / max_integral
                    #         ),
                    #     ],
                    #     [
                    #         1.0,
                    #         get_continuous_color(
                    #             COLORSCALE, intermed=max(0, n_photons[opid]) / max_integral
                    #         ),
                    #     ],
                    # ]

                    if ix % 2 == 0:
                        flip = 0
                    else:
                        flip = -1

                    opid_str = f"opid_{opid}"
                    light_plane = dict(
                        type="surface", 
                        y=np.full(xx.shape, det_bounds[ix][0][0][iside + flip])/10 - 240,
                        x=xx/10,
                        z=zz/10+1300,
                        opacity=0.4,
                        hoverinfo="text",
                        ids=[[opid_str, opid_str], [opid_str, opid_str]],
                        customdata=[[opid_str, opid_str], [opid_str, opid_str]],
                        # text=f"Optical detector {opid} waveform integral<br>{n_photons[opid]:.2e}",
                        # colorscale=light_color,
                        showlegend=False,
                        showscale=False,
                    )
                    drawn_objects.append(light_plane)
        return drawn_objects

    def construct_waveforms(self):
        waveforms = go.Figure()
        waveforms.update_xaxes(title_text='Time [ticks] (1 ns)')
        waveforms.update_yaxes(title_text='Adc counts')
        waveforms.update_layout(title_text='Waveform for optical detector')
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
                    # This div contains the plot on the left with 50% width of the viewport width
                    html.Div(dcc.Graph(
                        id='charge_light_tpc_plot',
                        style={'height': '80vh', 'width': '45vw'},  # Adjust size as needed
                        figure=self.tpc
                    ), style={'width': '45vw'}),  # This ensures the left plot takes up half the viewport width

                    # This div is the container for the two plots on the right
                    html.Div([
                        # This div contains the top right plot with 35% viewport width and 50% viewport height
                        html.Div(dcc.Graph(
                            id="light-waveform",
                            style={'height': '50vh', 'width': '35vw'},  # Adjust size as needed
                            figure=self.waveforms
                        )),
                        # This div contains the bottom right plot with 35% viewport width and 30% viewport height
                        html.Div(dcc.Graph(
                            id="another-graph",
                            style={'height': '30vh', 'width': '35vw'},  # Adjust size as needed
                            figure=self.larpix
                        )),
                    ], style={'display': 'flex', 'flexDirection': 'column', 'width': '45vw'}),
                ], style={'display': 'flex'})
            ]
        )
