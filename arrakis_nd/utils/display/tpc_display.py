import plotly.graph_objects as go
from dash import dcc,html
import numpy as np
import pandas as pd # remove this dependency later
import plotly
import h5py
import traceback

from arrakis_nd.utils.display.utils import custom_plotly_layout, get_continuous_color


class TPCDisplay:
    '''
    This class is used to visualize the events in three different plots:

    '''
    def __init__(
        self,
        id_suffix,
        geometry_info: dict = {}
    ):
        self.id_suffix = id_suffix
        self.interactions = None
        self.segments = None
        self.stack = None
        self.trajectories = None
        self.charge = None
        self.topology = None
        self.physics = None
        self.particle = None
        self.unique_topology = None
        self.geometry_info = {}
        self.hit_type = 'prompt'

        self.scale = 0.0
        self.plottype = 'q'
        self.topology_labels = {
            0: 'track',
            1: 'shower',
            2: 'blip'
        }
        self.topology_colors = {
            0: 'blue',
            1: 'green',
            2: 'orange',
        }
        self.physics_labels = {
            0: 'mip',
            1: 'hip',
            2: 'e-ionization',
            3: 'delta',
            4: 'michel',
            5: 'compton',
            6: 'conversion',
            7: 'nr',
            8: 'er'
        }
        self.physics_colors = {
            0: 'orange',
            1: 'magenta',
            2: 'blue',
            3: 'green',
            4: 'purple',
            5: 'yellow',
            6: 'maroon',
            7: 'deepskyblue',
            8: 'grey'
        }

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
                # self.light_traps = self.construct_light_detectors()
                # self.tpc.add_traces(self.light_traps)
            except Exception:
                print(traceback.format_exc())
                # print(f"issue with light traps {e}")

    def construct_detector(self):
        self.tpc = self.construct_tpcs()

    def construct_tpcs(self):
        tpc = go.Figure()
        tpc.update_layout(
            scene=dict(
                xaxis_title="x [cm]",
                yaxis_title="z [cm]",
                zaxis_title="y [cm]"
            ),
            autosize=True,
            margin=dict(l=0, r=0, b=0),  # Adjust margins to center the plot,
            legend=dict(
                x=0,
                y=1,
                xanchor='left',
                yanchor='top'
            ),
        )
        self.tpc_center, self.anodes, self.cathodes = self.draw_tpc()
        tpc.add_traces(self.tpc_center)
        tpc.add_traces(self.anodes)
        tpc.add_traces(self.cathodes)
        return tpc

    def construct_light_detectors(self, waveforms):
        """Plot optical detectors"""
        drawn_objects = []
        det_bounds = self.geometry_info["det_bounds"]
        channel_map = np.array([
            0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 65, 73, 81, 89, 97, 105, 113, 121,
            1, 9, 17, 25, 33, 41, 49, 57, 2, 10, 18, 26, 34, 42, 50, 58, 66, 74, 82, 90, 98, 106, 114, 122, 67,
            75, 83, 91, 99, 107, 115, 123, 3, 11, 19, 27, 35, 43, 51, 59, 4, 12, 20, 28, 36, 44, 52, 60, 68, 76,
            84, 92, 100, 108, 116, 124, 69, 77, 85, 93, 101, 109, 117, 125, 5, 13, 21, 29, 37, 45, 53, 61, 6, 14,
            22, 30, 38, 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 71, 79, 87, 95, 103, 111, 119, 127, 7, 15,
            23, 31, 39, 47, 55, 63
        ])  # this maps detector position to detector number
        # we need to invert the mapping because I'm stupid
        channel_map = np.argsort(channel_map)
        channel_map_deluxe = pd.read_csv('arrakis_nd/utils/display/sipm_channel_map.csv', header=0)

        xs = []
        ys = []
        zs = []
        for i in range(len(det_bounds)):
            if det_bounds[i][1] is True:
                xs.append([det_bounds[i][0][0][0], det_bounds[i][0][1][0]])
                ys.append([det_bounds[i][0][0][1], det_bounds[i][0][1][1]])
                zs.append([det_bounds[i][0][0][2], det_bounds[i][0][1][2]])

        COLORSCALE = plotly.colors.make_colorscale(
            plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.YlOrRd)[0]
        )

        photon_sums = []
        for i in range(len(xs)):
            opid = channel_map[i]
            # get all adc, channel belonging to opid=det_id.
            # We have a numpy array channel_map_deluxe with dtype
            # (det_id, tpc, side, sipm_pos, adc, channel)
            sipms = channel_map_deluxe[channel_map_deluxe['det_id'] == opid][['adc', 'channel']].values
            sum_photons = 0
            for adc, channel in sipms:
                wvfm = waveforms[:, int(adc), int(channel), :]
                sum_wvfm = np.sum(wvfm, axis=0)  # sum over the events
                sum_photons += np.sum(sum_wvfm, axis=0)  # sum over the time
            photon_sums.append(sum_photons)
        max_integral = np.max(photon_sums)

        for i in range(len(xs)):
            opid = channel_map[i]
            # get all adc, channel belonging to opid=det_id.
            # We have a numpy array channel_map_deluxe with dtype
            # (det_id, tpc, side, sipm_pos, adc, channel)
            sipms = channel_map_deluxe[channel_map_deluxe['det_id'] == opid][['adc', 'channel']].values
            sum_photons = 0
            for adc, channel in sipms:
                wvfm = waveforms[:, int(adc), int(channel), :]
                sum_wvfm = np.sum(wvfm, axis=0)
                sum_photons += np.sum(sum_wvfm, axis=0)
            opid_str = f"opid_{opid}"
            light_color = [
                [
                    0.0,
                    get_continuous_color(
                        COLORSCALE, intermed=max(0, sum_photons/max_integral)
                    ),
                ],
                [
                    1.0,
                    get_continuous_color(
                        COLORSCALE, intermed=max(0, sum_photons/max_integral)
                    ),
                ],
            ]
            light_plane = go.Surface(
                x=xs[i],
                z=ys[i],
                y=[zs[i], zs[i]],  # why flip y and z?
                colorscale=light_color,
                showscale=False,
                showlegend=False,
                opacity=0.2,
                hoverinfo="text",
                ids=[[opid_str, opid_str], [opid_str, opid_str]],
                text=f"Optical detector {opid} waveform integral<br>{max(0, sum_photons):.2e}",
            )

            drawn_objects.append(light_plane)
        self.tpc.add_traces(drawn_objects)

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
            z=[detector_center[1]],
            y=[detector_center[2]],
            marker=dict(size=3, color="green", opacity=0.5),
            mode="markers",
            name="tpc center",
        )
        anodes = self.draw_anode_planes(
            anode_xs,
            anode_ys,
            anode_zs,
            colorscale="ice",
            showscale=False,
            opacity=0.1
        )
        cathodes = self.draw_cathode_planes(
            anode_xs,
            anode_ys,
            anode_zs,
            colorscale="burg",
            showscale=False, opacity=0.1
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
        legendgroup = 'cathode_planes'  # Group name for legend
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
                showlegend = (i_z == 0 and i_x == 0)  # Show legend for the first trace only
                traces.append(go.Surface(
                    x=x,
                    y=z,
                    z=y,
                    legendgroup=legendgroup,
                    showlegend=showlegend,
                    name="Cathode Planes",
                    **kwargs
                ))

        return traces

    def draw_anode_planes(
        self,
        x_boundaries,
        y_boundaries,
        z_boundaries,
        **kwargs
    ):
        traces = []
        legendgroup = 'anode_planes'  # Group name for legend
        for i_z in range(int(len(z_boundaries) / 2)):
            for i_x in range(int(len(x_boundaries))):
                z, y = np.meshgrid(
                    np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2),
                    np.linspace(y_boundaries.min(), y_boundaries.max(), 2),
                )
                x = x_boundaries[i_x] * np.ones(z.shape)

                showlegend = (i_z == 0 and i_x == 0)  # Show legend for the first trace only
                traces.append(go.Surface(
                    x=x,
                    y=z,
                    z=y,
                    legendgroup=legendgroup,
                    showlegend=showlegend,
                    name="Anode Planes",
                    **kwargs
                ))

        return traces

    def update_flow_event(
        self,
        interactions,
        segments,
        stack,
        trajectories,
        charge,
    ):
        self.interactions = interactions
        self.segments = segments
        self.stack = stack
        self.trajectories = trajectories
        self.charge = charge

    def update_arrakis_event(
        self,
        topology,
        physics,
        particle,
        unique_topology
    ):
        self.topology = topology
        self.physics = physics
        self.particle = particle
        self.unique_topology = unique_topology

    def update_plottype(
        self,
        plottype
    ):
        self.plottype = plottype

    def plot_event(self):
        if self.charge is not None:
            if self.scale == 0.0:
                marker_size = 5.0  # Fixed size when scale is 0
            else:
                marker_size = np.array([q * self.scale for q in self.charge['Q']])  # Scale marker size

            names_to_remove = [
                f'calib_{self.hit_type}_hits',
                'track', 'shower', 'blip',
                'mip', 'hip', 'e-ionization', 'delta', 'michel', 'compton', 'conversion', 'nr', 'er'
            ]
            # Create a new list of traces that excludes the ones with the specified name
            new_traces = [trace for trace in self.tpc.data if trace.name not in names_to_remove]
            # Update the figure with the new list of traces
            self.tpc.data = new_traces

            if self.plottype == 'q':
                traces = self.plot_q(marker_size)
            elif self.plottype == 'topology':
                traces = self.plot_topology(marker_size)
            elif self.plottype == 'physics':
                traces = self.plot_physics(marker_size)

            self.tpc.add_traces(traces)

    def highlight_point(self, hit_id, color='red'):
        """
        Highlight the specified point by setting its color and opacity.
        """
        # Default color and opacity for all points
        default_opacity = 0.5
        highlighted_opacity = 1.0

        # Update traces to highlight the selected point
        new_data = []
        for trace in self.tpc.data:
            customdata = trace['customdata']

            if 'customdata' in trace and customdata is not None:
                updated_marker = trace['marker']
                updated_marker['opacity'] = [
                    highlighted_opacity if customdata[i][0] == hit_id else default_opacity
                    for i in range(len(customdata))
                ]
                updated_marker['color'] = [
                    color if customdata[i][0] == hit_id else updated_marker['color']
                    for i in range(len(customdata))
                ]

                new_trace = trace.update(marker=updated_marker)
                new_data.append(new_trace)
            else:
                new_data.append(trace)

        self.tpc.data = new_data

    def plot_q(
        self,
        marker_size
    ):
        charge_hits_traces = go.Scatter3d(
            x=self.charge['x'],
            z=self.charge['y'],
            y=self.charge['z'],
            marker={
                "size": marker_size,
                "opacity": 0.5,
                "colorscale": "cividis",
            },
            name=f"calib_{self.hit_type}_hits",
            mode="markers",
            showlegend=True,
            opacity=0.5,
            customdata=np.array(list(zip(
                list(range(1, len(self.charge['x']) + 1)),
                self.charge['Q'],
                self.charge['E'],
                self.particle,
                self.unique_topology
            ))),
            hovertemplate=(
                '<b>hit_id:</b> %{customdata[0]}<br>'
                '<b>x:</b> %{x:.3f}<br>'
                '<b>y:</b> %{y:.3f}<br>'
                '<b>z:</b> %{z:.3f}<br>'
                '<b>Q:</b> %{customdata[1]:.3f}<br>'
                '<b>E:</b> %{customdata[2]:.3f}<br>'
                '<b>pdg_id:</b> %{customdata[3]}<br>'
                '<b>track_id:</b> %{customdata[4]}<extra></extra>'
            )
        )
        return charge_hits_traces

    def plot_topology(
        self,
        marker_size
    ):
        traces = []
        for unique_topology in np.unique(self.topology):
            topology_mask = (self.topology == unique_topology)
            traces.append(go.Scatter3d(
                x=self.charge['x'][topology_mask],
                z=self.charge['y'][topology_mask],
                y=self.charge['z'][topology_mask],
                marker={
                    "size": marker_size,
                    "opacity": 0.5,
                    "color": self.topology_colors[unique_topology]
                },
                name=self.topology_labels[unique_topology],
                mode="markers",
                showlegend=True,
                opacity=0.5,
                customdata=np.array(list(zip(
                    np.array(list(range(1, len(self.charge['x']) + 1)))[topology_mask],
                    self.charge['Q'][topology_mask],
                    self.charge['E'][topology_mask],
                    self.particle[topology_mask],
                    self.unique_topology[topology_mask]
                ))),
                hovertemplate=(
                    '<b>hit_id:</b> %{customdata[0]}<br>'
                    '<b>x:</b> %{x:.3f}<br>'
                    '<b>y:</b> %{y:.3f}<br>'
                    '<b>z:</b> %{z:.3f}<br>'
                    '<b>Q:</b> %{customdata[1]:.3f}<br>'
                    '<b>E:</b> %{customdata[2]:.3f}<br>'
                    '<b>pdg_id:</b> %{customdata[3]}<br>'
                    '<b>track_id:</b> %{customdata[4]}<extra></extra>'
                )
            ))
        return traces

    def plot_physics(
        self,
        marker_size
    ):
        traces = []
        for unique_physics in np.unique(self.physics):
            physics_mask = (self.physics == unique_physics)
            traces.append(go.Scatter3d(
                x=self.charge['x'][physics_mask],
                z=self.charge['y'][physics_mask],
                y=self.charge['z'][physics_mask],
                marker={
                    "size": marker_size,
                    "opacity": 0.5,
                    "color": self.physics_colors[unique_physics]
                },
                name=self.physics_labels[unique_physics],
                mode="markers",
                showlegend=True,
                opacity=0.5,
                customdata=np.array(list(zip(
                    np.array(list(range(1, len(self.charge['x']) + 1)))[physics_mask],
                    self.charge['Q'][physics_mask],
                    self.charge['E'][physics_mask],
                    self.particle[physics_mask],
                    self.unique_topology[physics_mask]
                ))),
                hovertemplate=(
                    '<b>hit_id:</b> %{customdata[0]}<br>'
                    '<b>x:</b> %{x:.3f}<br>'
                    '<b>y:</b> %{y:.3f}<br>'
                    '<b>z:</b> %{z:.3f}<br>'
                    '<b>Q:</b> %{customdata[1]:.3f}<br>'
                    '<b>E:</b> %{customdata[2]:.3f}<br>'
                    '<b>pdg_id:</b> %{customdata[3]}<br>'
                    '<b>track_id:</b> %{customdata[4]}<extra></extra>'
                )
            ))
        return traces

    def generate_layout(self):
        self.construct_detector()
        self.plot_event()
        self.layout = html.Div([
            dcc.Graph(
                id=f'tpc_plot_{self.id_suffix}',
                style={'height': '60vh', 'width': '40vw'},
                figure=self.tpc,
                config={'displayModeBar': True}
            )],
            style={'width': '40vw'}
        )
