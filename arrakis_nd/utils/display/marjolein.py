# In my_dash_app.py
from dash import Dash
import dash_core_components as dcc
import dash_uploader as du
import dash_html_components as html
import plotly.graph_objects as go
import atexit
import plotly
import plotly.graph_objects as go
import shutil

from dash import dcc
from dash import html
#from dash import no_update
from dash.exceptions import PreventUpdate
from dash_extensions.enrich import Output, DashProxy, Input, State, MultiplexerTransform

# from display_utils import parse_contents, create_3d_figure, plot_waveform

from os.path import basename
from pathlib import Path










# # Settings and constants
# UPLOAD_FOLDER_ROOT = "cache"
# from dash.dependencies import Input, Output
# # plus other necessary imports

# def parse_contents(filename):
#     data = h5flow.data.H5FlowDataManager(filename, "r")
#     num_events = data["charge/events/data"].shape[0]
#     return data, num_events


# def create_3d_figure(data, evid):
#     fig = go.Figure()
#     print("here we go")
#     # Select the hits for the current event
#     prompthits_ev = data["charge/events", "charge/calib_prompt_hits", evid]
#     finalhits_ev = data["charge/events", "charge/calib_final_hits", evid]
#     # select the segments (truth) for the current event
#     try:
#         prompthits_segs = data[
#             "charge/events",
#             "charge/calib_prompt_hits",
#             "charge/packets",
#             "mc_truth/segments",  # called segments in minirun4
#             evid,
#         ]
#         sim_version = "minirun4"
#         print("Found truth info in minirun4 format")
#     except:
#         print("No truth info in minirun4 format found")
#         try:
#             prompthits_segs = data[
#                 "charge/events",
#                 "charge/calib_prompt_hits",
#                 "charge/packets",
#                 "mc_truth/tracks",  # called tracks in minirun3
#                 evid,
#             ]
#             sim_version = "minirun3"
#             print("Found truth info in minirun3 format")
#         except:
#             print("No truth info in minirun3 format found")
#             prompthits_segs = None

#     # Plot the prompt hits
#     print("Plotting prompt hits")
#     prompthits_traces = go.Scatter3d(
#         x=prompthits_ev.data["x"].flatten(),
#         y=(prompthits_ev.data["y"].flatten()),
#         z=(prompthits_ev.data["z"].flatten()),
#         marker_color=prompthits_ev.data["E"].flatten()
#         * 1000,  # convert to MeV from GeV for minirun4, not sure for minirun3
#         marker={
#             "size": 1.75,
#             "opacity": 0.7,
#             "colorscale": "cividis",
#             "colorbar": {
#                 "title": "Hit energy [MeV]",
#                 "titlefont": {"size": 12},
#                 "tickfont": {"size": 10},
#                 "thickness": 15,
#                 "len": 0.5,
#                 "xanchor": "left",
#                 "x": 0,
#             },
#         },
#         name="prompt hits",
#         mode="markers",
#         showlegend=True,
#         opacity=0.7,
#         customdata=prompthits_ev.data["E"].flatten() * 1000,
#         hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
#     )
#     print("Adding prompt hits to figure")
#     fig.add_traces(prompthits_traces)

#     # Plot the final hits
#     print("Plotting final hits")
#     finalhits_traces = go.Scatter3d(
#         x=finalhits_ev.data["x"].flatten(),
#         y=(finalhits_ev.data["y"].flatten()),
#         z=(finalhits_ev.data["z"].flatten()),
#         marker_color=finalhits_ev.data["E"].flatten() * 1000,
#         marker={
#             "size": 1.75,
#             "opacity": 0.7,
#             "colorscale": "Plasma",
#             "colorbar": {
#                 "title": "Hit energy [MeV]",
#                 "titlefont": {"size": 12},
#                 "tickfont": {"size": 10},
#                 "thickness": 15,
#                 "len": 0.5,
#                 "xanchor": "left",
#                 "x": 0,
#             },
#         },
#         name="final hits",
#         mode="markers",
#         visible="legendonly",
#         showlegend=True,
#         opacity=0.7,
#         customdata=finalhits_ev.data["E"].flatten() * 1000,
#         hovertemplate="<b>x:%{x:.3f}</b><br>y:%{y:.3f}<br>z:%{z:.3f}<br>E:%{customdata:.3f}",
#     )
#     print("Adding final hits to figure")
#     fig.add_traces(finalhits_traces)

#     print("Plotting segments")
#     if prompthits_segs is not None:
#         segs_traces = plot_segs(
#             prompthits_segs[0, :, 0, 0],
#             sim_version=sim_version,
#             mode="lines",
#             name="edep segments",
#             visible="legendonly",
#             line_color="red",
#             showlegend=True,
#         )
#         print("Adding segments to figure")
#         fig.add_traces(segs_traces)

#     # Draw the TPC
#     print("Drawing TPC")
#     tpc_center, anodes, cathodes = draw_tpc(sim_version)
#     light_detectors = draw_light_detectors(data, evid)

#     print("Adding TPC to figure")
#     fig.add_traces(tpc_center)
#     fig.add_traces(anodes)
#     fig.add_traces(cathodes)
#     print("Adding light detectors to figure")
#     fig.add_traces(light_detectors)

#     return fig


# def plot_segs(segs, sim_version="minirun4", **kwargs):
#     def to_list(axis):
#         if sim_version == "minirun4":
#             nice_array = np.column_stack(
#                 [segs[f"{axis}_start"], segs[f"{axis}_end"], np.full(len(segs), None)]
#             ).flatten()
#         if sim_version == "minirun3":
#             nice_array = np.column_stack(
#                 [
#                     segs[f"{axis}_start"] * 10,
#                     segs[f"{axis}_end"] * 10,
#                     np.full(len(segs), None),
#                 ]
#             ).flatten()
#         return nice_array

#     x, y, z = (to_list(axis) for axis in "xyz")

#     trace = go.Scatter3d(x=x, y=y, z=z, **kwargs)

#     return trace


# def draw_tpc(sim_version="minirun4"):
#     anode_xs = np.array([-63.931, -3.069, 3.069, 63.931])
#     anode_ys = np.array([-19.8543, 103.8543])  # two ys
#     anode_zs = np.array([-64.3163, -2.6837, 2.6837, 64.3163])  # four zs
#     if sim_version == "minirun4":  # hit coordinates are in cm
#         detector_center = (0, -268, 1300)
#         anode_ys = anode_ys - (268 + 42)
#         anode_zs = anode_zs + 1300
#     if sim_version == "minirun3":  # hit coordinates are in mm
#         detector_center = (0, 42 * 10, 0)
#         anode_xs = anode_xs * 10
#         anode_ys = anode_ys * 10
#         anode_zs = anode_zs * 10

#     center = go.Scatter3d(
#         x=[detector_center[0]],
#         y=[detector_center[1]],
#         z=[detector_center[2]],
#         marker=dict(size=3, color="green", opacity=0.5),
#         mode="markers",
#         name="tpc center",
#     )
#     anodes = draw_anode_planes(
#         anode_xs, anode_ys, anode_zs, colorscale="ice", showscale=False, opacity=0.1
#     )
#     cathodes = draw_cathode_planes(
#         anode_xs, anode_ys, anode_zs, colorscale="burg", showscale=False, opacity=0.1
#     )
#     return center, anodes, cathodes


# def draw_cathode_planes(x_boundaries, y_boundaries, z_boundaries, **kwargs):
#     traces = []
#     for i_z in range(int(len(z_boundaries) / 2)):
#         for i_x in range(int(len(x_boundaries) / 2)):
#             z, y = np.meshgrid(
#                 np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2),
#                 np.linspace(y_boundaries.min(), y_boundaries.max(), 2),
#             )
#             x = (
#                 (x_boundaries[i_x * 2] + x_boundaries[i_x * 2 + 1])
#                 * 0.5
#                 * np.ones(z.shape)
#             )
#             traces.append(go.Surface(x=x, y=y, z=z, **kwargs))

#     return traces


# def draw_anode_planes(x_boundaries, y_boundaries, z_boundaries, **kwargs):
#     traces = []
#     for i_z in range(int(len(z_boundaries) / 2)):
#         for i_x in range(int(len(x_boundaries))):
#             z, y = np.meshgrid(
#                 np.linspace(z_boundaries[i_z * 2], z_boundaries[i_z * 2 + 1], 2),
#                 np.linspace(y_boundaries.min(), y_boundaries.max(), 2),
#             )
#             x = x_boundaries[i_x] * np.ones(z.shape)

#             traces.append(go.Surface(x=x, y=y, z=z, **kwargs))

#     return traces


# def draw_light_detectors(data, evid):
#     try:
#         charge = data["charge/events", evid][["id", "unix_ts"]]
#         num_light = data["light/events/data"].shape[0]
#         light = data["light/events", slice(0, num_light)][
#             ["id", "utime_ms"]
#         ]  # we have to try them all, events may not be time ordered
#     except:
#         print("No light information found, not plotting light detectors")
#         return []

#     match_light = match_light_to_charge_event(charge, light, evid)

#     if match_light is None:
#         print(
#             f"No light event matches found for charge event {evid}, not plotting light detectors"
#         )
#         return []

#     waveforms_all_detectors = get_waveforms_all_detectors(data, match_light)

#     # make a list of the sum of the waveform and the channel index
#     integral = np.sum(np.sum(waveforms_all_detectors, axis=2), axis=0)
#     max_integral = np.max(integral)
#     index = np.arange(0, waveforms_all_detectors.shape[1], 1)

#     # plot for each of the 96 channels per tpc the sum of the adc values
#     drawn_objects = []
#     drawn_objects.extend(plot_light_traps(data, integral, index, max_integral))

#     return drawn_objects


# def match_light_to_charge_event(charge, light, evid):
#     """
#     Match the light events to the charge event by looking at proximity in time.
#     Use unix time for this, since it should refer to the same time in both readout systems.
#     For now we just take all the light within 1s from the charge event time.
#     """
#     matches = []
#     for i in range(len(light)):
#         if np.abs(light["utime_ms"][i][0] / 1000 - charge["unix_ts"][0]) < 0.5:
#             matches.append([charge["id"][0], light["id"][i]])

#     match_light = []
#     for i in range(len(matches)):
#         if (
#             matches[i][0] == evid
#         ):  # just checking that we get light for the right charge event
#             match_light.append(matches[i][1])
#     if len(match_light) == 0:
#         match_light = None  # no light for this charge event

#     return match_light


# def get_waveforms_all_detectors(data, match_light):
#     """
#     Get the light waveforms for the matched light events.
#     """
#     light_wvfm = data["/light/wvfm", match_light]

#     samples_mod0 = light_wvfm["samples"][:, 0:2, :, :]
#     samples_mod1 = light_wvfm["samples"][:, 2:4, :, :]
#     samples_mod2 = light_wvfm["samples"][:, 4:6, :, :]
#     samples_mod3 = light_wvfm["samples"][:, 6:8, :, :]

#     sipm_channels_module0 = np.array(
#         [2, 3, 4, 5, 6, 7]
#         + [9, 10]
#         + [11, 12]
#         + [13, 14]
#         + [18, 19, 20, 21, 22, 23]
#         + [25, 26]
#         + [27, 28]
#         + [29, 30]
#         + [34, 35, 36, 37, 38, 39]
#         + [41, 42]
#         + [43, 44]
#         + [45, 46]
#         + [50, 51, 52, 53, 54, 55]
#         + [57, 58]
#         + [59, 60]
#         + [61, 62]
#     )

#     sipm_channels_modules = np.array(
#         [4, 5, 6, 7, 8, 9]
#         + [10, 11, 12, 13, 14, 15]
#         + [20, 21, 22, 23, 24, 25]
#         + [26, 27, 28, 29, 30, 31]
#         + [36, 37, 38, 39, 40, 41]
#         + [42, 43, 44, 45, 46, 47]
#         + [52, 53, 54, 55, 56, 57]
#         + [58, 59, 60, 61, 62, 63]
#     )
#     adcs_mod0 = samples_mod0[:, :, sipm_channels_module0, :]
#     adcs_mod1 = samples_mod1[:, :, sipm_channels_modules, :]
#     adcs_mod2 = samples_mod2[:, :, sipm_channels_modules, :]
#     adcs_mod3 = samples_mod3[:, :, sipm_channels_modules, :]

#     all_adcs = np.concatenate((adcs_mod0, adcs_mod1, adcs_mod2, adcs_mod3), axis=1)

#     # instead of a (m, 8, 48, 1000) array, we want a (m, 4, 96, 1000) array
#     # modules instead of tpcs, and 96 channels per module
#     m = len(match_light)
#     all_modules = all_adcs.reshape((m, 4, 96, 1000))

#     # now we make a full array for all the modules
#     # could have been done in one step, but this is easier to read
#     all_detector = all_modules.reshape((m, 384, 1000))

#     return all_detector


# def plot_light_traps(data, n_photons, op_indeces, max_integral):
#     """Plot optical detectors"""
#     drawn_objects = []
#     ys = np.flip(
#         np.array(
#             [
#                 -595.43,
#                 -545.68,
#                 -490.48,
#                 -440.73,
#                 -385.53,
#                 -335.78,
#                 -283.65,
#                 -236.65,
#                 -178.70,
#                 -131.70,
#                 -73.75,
#                 -26.75,
#                 25.38,
#                 75.13,
#                 130.33,
#                 180.08,
#                 235.28,
#                 285.03,
#                 337.15,
#                 384.15,
#                 442.10,
#                 489.10,
#                 547.05,
#                 594.05,
#             ]
#         )
#         / 10
#     )
#     light_width = ys[1] - ys[0]

#     det_bounds = data["/geometry_info/det_bounds/data"]
#     COLORSCALE = plotly.colors.make_colorscale(
#         plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.YlOrRd)[0]
#     )

#     for ix in range(0, det_bounds.shape[0]):
#         for ilight, light_y in enumerate(ys):
#             for iside in range(2):
#                 opid = ilight + iside * len(ys) + ix * len(ys) * 2
#                 if opid not in op_indeces:
#                     continue
#                 xx = np.linspace(det_bounds[ix][0][0][0], det_bounds[ix][0][1][0], 2)
#                 zz = np.linspace(
#                     light_y - light_width / 2 + det_bounds[0][0][1] + 0.25,
#                     light_y + light_width / 2 + det_bounds[0][0][1] - 0.25,
#                     2,
#                 )

#                 xx, zz = np.meshgrid(xx, zz)
#                 light_color = [
#                     [
#                         0.0,
#                         get_continuous_color(
#                             COLORSCALE, intermed=max(0, n_photons[opid]) / max_integral
#                         ),
#                     ],
#                     [
#                         1.0,
#                         get_continuous_color(
#                             COLORSCALE, intermed=max(0, n_photons[opid]) / max_integral
#                         ),
#                     ],
#                 ]

#                 if ix % 2 == 0:
#                     flip = 0
#                 else:
#                     flip = -1

#                 opid_str = f"opid_{opid}"
#                 light_plane = dict(
#                     type="surface", 
#                     y=np.full(xx.shape, det_bounds[ix][0][0][iside + flip])/10 -240,
#                     x=xx/10,
#                     z=zz/10+1300,
#                     opacity=0.4,
#                     hoverinfo="text",
#                     ids=[[opid_str, opid_str], [opid_str, opid_str]],
#                     customdata=[[opid_str, opid_str], [opid_str, opid_str]],
#                     text=f"Optical detector {opid} waveform integral<br>{n_photons[opid]:.2e}",
#                     colorscale=light_color,
#                     showlegend=False,
#                     showscale=False,
#                 )

#                 drawn_objects.append(light_plane)

#     return drawn_objects

# def plot_waveform(data, evid, opid):
#     try:
#         charge = data["charge/events", evid][["id", "unix_ts"]]
#         num_light = data["light/events/data"].shape[0]
#         light = data["light/events", slice(0, num_light)][
#             ["id", "utime_ms"]
#         ]  # we have to try them all, events may not be time ordered
#     except:
#         print("No light information found, not plotting light waveform")
#         return []

#     match_light = match_light_to_charge_event(charge, light, evid)

#     if match_light is None:
#         print(
#             f"No light event matches found for charge event {evid}, not plotting light waveform"
#         )
#         return []

#     fig = go.Figure()
#     waveforms_all_detectors = get_waveforms_all_detectors(data, match_light)
#     wvfm_opid = waveforms_all_detectors[:, opid, :]
    
#     x = np.arange(0, 1000, 1)
#     y = np.sum(wvfm_opid, axis=0)
#     drawn_objects = go.Scatter(x=x, y=y)
#     fig.add_traces(drawn_objects)
    
#     fig.update_xaxes(title_text='Time [ticks] (1 ns)')
#     fig.update_yaxes(title_text='Adc counts')
#     fig.update_layout(title_text=f'Waveform for optical detector {opid}')
#     return fig



# def get_continuous_color(colorscale, intermed):
#     """
#     Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
#     color for any value in that range.

#     Plotly doesn't make the colorscales directly accessible in a common format.
#     Some are ready to use:

#         colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

#     Others are just swatches that need to be constructed into a colorscale:

#         viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
#         colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

#     :param colorscale: A plotly continuous colorscale defined with RGB string colors.
#     :param intermed: value in the range [0, 1]
#     :return: color in rgb string format
#     :rtype: str
#     """
#     if len(colorscale) < 1:
#         raise ValueError("colorscale must have at least one color")

#     if intermed <= 0 or len(colorscale) == 1:
#         return colorscale[0][1]
#     if intermed >= 1:
#         return colorscale[-1][1]

#     for cutoff, color in colorscale:
#         if intermed > cutoff:
#             low_cutoff, low_color = cutoff, color
#         if intermed <= cutoff:
#             high_cutoff, high_color = cutoff, color
#             break

#     # noinspection PyUnboundLocalVariable
#     return plotly.colors.find_intermediate_color(
#         lowcolor=low_color,
#         highcolor=high_color,
#         intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
#         colortype="rgb",
#     )

# def run_my_dash_app():
#     # Create the app
#     app = DashProxy(__name__, title="2x2 event display")
#     du.configure_upload(app, UPLOAD_FOLDER_ROOT)  # without this upload will not work

#     # App layout
#     app.layout = html.Div(
#         [
#             # Hidden divs to store data
#             dcc.Location(id="url"),
#             dcc.Store(id="filename", storage_type="local", data=None),
#             dcc.Store(id='data-length', data=0),
#             # Header
#             html.H1(children="2x2 event display", style={"textAlign": "center"}),
#             html.Div(children="", id="filename-div", style={"textAlign": "center"}),
#             # Upload button
#             html.Div(
#                 du.Upload(
#                     id="upload-data-div",
#                     text="Upload Flow HDF5 File",
#                     max_file_size=10000,
#                     chunk_size=5,
#                     default_style={
#                         "width": "15em",
#                         "padding": "0",
#                         "margin": "0",
#                     },
#                     pause_button=True,
#                     filetypes=["h5"],
#                 ),
#             ),
#             # Event ID input box
#             dcc.Input(
#                     id="input-evid",
#                     type="number",
#                     placeholder="0",
#                     debounce=True,
#                     style={
#                         "width": "6em",
#                         "display": "inline-block",
#                         "margin-right": "0.5em",
#                         "margin-left": "0.5em",
#                     },
#                 ),
#             # Event ID buttons
#             html.Button('Previous Event', id='prev-button', n_clicks=0),
#             html.Button('Next Event', id='next-button', n_clicks=0),
#             dcc.Store(id='event-id', data=0),
#             html.Div(id='evid-div', style={"textAlign": "center"}),
#             # Graphs
#             html.Div([
#                 # Existing 3D graph
#                 html.Div(dcc.Graph(id='3d-graph', style={'height': '70vh', 'width': '50vw'})),

#                 # New Light waveform graph
#                 html.Div(dcc.Graph(id="light-waveform", style={'height': '50vh', 'width': '35vw'})),

                
#             ], style={'display': 'flex'}),
#             # New Another Graph (replace with your actual component)
#             html.Div(dcc.Graph(id="another-graph", style={'height': '30vh', 'width': '35vw', 'float': 'right'})),
#         ]
#     )


#     # Callbacks


#     # Callback to handle the upload
#     # =============================
#     @app.callback(
#         [
#             Output("filename", "data"),
#             Output("filename-div", "children"),
#             Output("event-id", "data", allow_duplicate=True),
#             Output('data-length', 'data'),
#         ],
#         [
#             Input("upload-data-div", "isCompleted"),
#         ],
#         [
#             State("filename", "data"),
#             State("upload-data-div", "fileNames"),
#             State("upload-data-div", "upload_id"),
#         ],
#         prevent_initial_call=True
#     )
#     def upload_file(is_completed, current_filename, filenames, upload_id):
#         """
#         Upload HDF5 file to cache. If the upload is completed,
#         update the filename. Initialise the event ID to 0.
#         """
#         if not is_completed:
#             raise PreventUpdate

#         if filenames is not None:
#             if upload_id:
#                 root_folder = Path(UPLOAD_FOLDER_ROOT) / upload_id
#             else:
#                 root_folder = Path(UPLOAD_FOLDER_ROOT)
#             _, num_events = parse_contents(str(root_folder / filenames[0]))
#             new_filename = str(root_folder / filenames[0])
#             return new_filename, basename(filenames[0]), 0, num_events

#         return current_filename, "no file uploaded", 0, 0



#     # Callbacks to handle the event ID
#     # =================================
#     @app.callback(
#         Output('event-id', 'data', allow_duplicate=True),
#         Input('next-button', 'n_clicks'),
#         State('event-id', 'data'),
#         State('data-length', 'data'),
#         prevent_initial_call=True
#     )
#     def increment(n, evid, max_value):
#         if n > 0:
#             new_evid = evid + 1
#             if new_evid > max_value:  # wrap around
#                 return 0
#             else:
#                 return new_evid

#     @app.callback(
#         Output('event-id', 'data', allow_duplicate=True),
#         Input('prev-button', 'n_clicks'),
#         State('event-id', 'data'),
#         State('data-length', 'data'),
#         prevent_initial_call=True
#     )
#     def decrement(n, evid, max_value):
#         if n > 0:
#             if evid > 0:
#                 return evid - 1
#             return max_value - 1 # wrap around

#     @app.callback(
#         Output('event-id', 'data', allow_duplicate=True),
#         Input('input-evid', 'value'),
#         State('data-length', 'data'),
#         prevent_initial_call=True
#     )
#     def set_evid(value, max_value):
#         if value is not None:
#             if value > max_value:  # not possible to go higher than max value
#                 return max_value
#             else:
#                 return value

#     @app.callback(
#         Output('evid-div', 'children'),
#         Input('event-id', 'data'),
#         State('data-length', 'data'),
#     )
#     def update_div(evid, max_value):
#         return f'Event ID: {evid}/{max_value}'

#     # Callback to display the event
#     # =============================
#     @app.callback(
#         Output('3d-graph', 'figure'),
#         Input('filename', 'data'),
#         Input('event-id', 'data'),
#         prevent_initial_call=True
#     )
#     def update_graph(filename, evid):
#         if filename is not None:
#             data, _ = parse_contents(filename)
#             return create_3d_figure(data, evid)
        
#     @app.callback(
#         Input('filename', 'data'),
#         Input('event-id', 'data'),
#         Input('3d-graph', 'figure'),
#         Input('3d-graph', 'clickData'),
#         Output('light-waveform', 'figure'),
#     )
#     def update_light_waveform(filename, evid, graph, click_data):
#         if click_data:
#             curvenum = int(click_data["points"][0]["curveNumber"])
#             print(curvenum) # this is the curve number from the clickdata of the graph
#             opid = int(graph['data'][curvenum]['ids'][0][0].split('_')[1])
#             print(opid) # the opid related to the curvenumber
#             if filename is not None:
#                 data, _ = parse_contents(filename)
#                 return plot_waveform(data, evid, opid)
#         return go.Figure()


#     # Cleaning up
#     # ===========
#     @atexit.register
#     def clean_cache():
#         """Delete uploaded files"""
#         try:
#             shutil.rmtree(Path(UPLOAD_FOLDER_ROOT))
#         except OSError as err:
#             print("Can't clean %s : %s" % (UPLOAD_FOLDER_ROOT, err.strerror))

#     app.run(jupyter_mode="jupyterlab")