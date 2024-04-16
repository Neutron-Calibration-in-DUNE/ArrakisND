import plotly
import plotly.graph_objects as go


def custom_legend_name(fig,new_names):
    '''
    Custom legend name for plotly figures

    Args:
        fig    (px.Figure): plotly figure
        new_names (list):      list of new names for the legend

    Returns:
        fig (px.Figure): plotly figure with the custom legend names
    '''
    for i, new_name in enumerate(new_names):
        fig.data[i].name = new_name
        # fig.data[i].hovertemplate = new_name
    # fig.for_each_trace(lambda t: t.update(name = dict[t.name],legendgroup = dict[t.name],hovertemplate = t.hovertemplate.replace(t.name, dict[t.name])))
    return fig

def custom_plotly_layout(fig,xaxis_title ="",yaxis_title="",title="",legend_pos=["out","h"],log_menu=True,fontsize=16,figsize=None,ranges=(None,None),matches=("x","y"),
    tickformat=('.s','.s'),barmode="stack",margin={"auto":True}):
    '''
    Custom layout for plotly figures

    Args:
        fig      (Figure):       plotly figure
        xaxis_title (str):       x axis title (default: "")
        yaxis_title (str):       y axis title (default: "")
        title       (str):       title (default: "")
        legend_pos  (list):      legend position (default: ["out","h"])
        log_menu    (bool):      log menu (default: True)
        fontsize    (int):       font size (default: 16)
        figsize     (tuple):     figure size (default: None)
        ranges      (tuple):     axis ranges (default: (None,None))
        matches     (tuple):     axis matches (default: ("x","y"))
        tickformat  (tuple):     axis tick format (default: ('.s','.s'))
        barmode     (str):       barmode (default: "stack")
        margin      (dict):      figure margin (default: {"auto":True,"color":"white","margin":(0,0,0,0)})
    
    Returns:
        fig (Figure): plotly figure with the custom layout
    '''
    default_margin = {"color":"white","margin":(0,0,0,0)}
    if margin != None:
            for key in default_margin.keys():
                if key not in margin.keys():
                    margin[key] = default_margin[key]

    legend_opt = {"in": dict(yanchor="top", xanchor="right", x=0.99),"out": dict(yanchor="bottom", xanchor="right", x=1.2,y=0.2) }
    if log_menu: fig.update_layout( updatemenus=[ dict( buttons=list([ dict(args=[{"xaxis.type": "linear", "yaxis.type": "linear"}], label="LinearXY", method="relayout"),
                                                             dict(args=[{"xaxis.type": "log", "yaxis.type": "log"}],       label="LogXY",    method="relayout"),
                                                             dict(args=[{"xaxis.type": "linear", "yaxis.type": "log"}],    label="LogY",     method="relayout"),
                                                             dict(args=[{"xaxis.type": "log", "yaxis.type": "linear"}],    label="LogX",     method="relayout") ]),
                          direction="down", pad={"r": 10, "t": 10}, showactive=True, x=-0.1, xanchor="left", y=1.5, yanchor="top" ) ] )
    fig.update_layout( template="presentation", title=title, xaxis_title=xaxis_title, yaxis_title=yaxis_title, barmode=barmode,
                          font=dict(size=fontsize,family="serif"), legend_title_text='', legend = legend_opt[legend_pos[0]],legend_orientation=legend_pos[1], showlegend=True, paper_bgcolor=margin["color"])
    fig.update_xaxes(matches=matches[0],showline=True,mirror="ticks",zeroline=False,showgrid=True,minor_ticks="inside",tickformat=tickformat[0],range=ranges[0])
    fig.update_yaxes(matches=matches[1],showline=True,mirror="ticks",zeroline=False,showgrid=True,minor_ticks="inside",tickformat=tickformat[1],range=ranges[1])

    if figsize != None: fig.update_layout(width=figsize[0],height=figsize[1])
    if margin["auto"] == False: fig.update_layout(margin=dict(l=margin["margin"][0], r=margin["margin"][1], t=margin["margin"][2], b=margin["margin"][3]))

    return fig

def get_continuous_color(colorscale, intermed):
    """
    Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    color for any value in that range.

    Plotly doesn't make the colorscales directly accessible in a common format.
    Some are ready to use:

        colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    Others are just swatches that need to be constructed into a colorscale:

        viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
        colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    :param intermed: value in the range [0, 1]
    :return: color in rgb string format
    :rtype: str
    """
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    if intermed <= 0 or len(colorscale) == 1:
        return colorscale[0][1]
    if intermed >= 1:
        return colorscale[-1][1]

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        if intermed <= cutoff:
            high_cutoff, high_color = cutoff, color
            break

    # noinspection PyUnboundLocalVariable
    return plotly.colors.find_intermediate_color(
        lowcolor=low_color,
        highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb",
    )

# OLD FUNCTIONS (TO BE REMOVED)
# def plot_segs(segs, sim_version="minirun4", **kwargs):
#         def to_list(axis):
#             if sim_version == "minirun4":
#                 nice_array = np.column_stack(
#                     [segs[f"{axis}_start"], segs[f"{axis}_end"], np.full(len(segs), None)]
#                 ).flatten()
#             if sim_version == "minirun3":
#                 nice_array = np.column_stack(
#                     [
#                         segs[f"{axis}_start"] * 10,
#                         segs[f"{axis}_end"] * 10,
#                         np.full(len(segs), None),
#                     ]
#                 ).flatten()
#             return nice_array

#         x, y, z = (to_list(axis) for axis in "xyz")

#         trace = go.Scatter3d(x=x, y=y, z=z, **kwargs)

#         return trace  
    
# def draw_tpc(sim_version="minirun4"):
#         anode_xs = np.array([-63.931, -3.069, 3.069, 63.931])
#         anode_ys = np.array([-19.8543, 103.8543])  # two ys
#         anode_zs = np.array([-64.3163, -2.6837, 2.6837, 64.3163])  # four zs
#         if sim_version == "minirun4":  # hit coordinates are in cm
#             detector_center = (0, -268, 1300)
#             anode_ys = anode_ys - (268 + 42)
#             anode_zs = anode_zs + 1300
#         if sim_version == "minirun3":  # hit coordinates are in mm
#             detector_center = (0, 42 * 10, 0)
#             anode_xs = anode_xs * 10
#             anode_ys = anode_ys * 10
#             anode_zs = anode_zs * 10

#         center = go.Scatter3d(
#             x=[detector_center[0]],
#             y=[detector_center[1]],
#             z=[detector_center[2]],
#             marker=dict(size=3, color="green", opacity=0.5),
#             mode="markers",
#             name="tpc center",
#         )
#         anodes = draw_anode_planes(
#             anode_xs, anode_ys, anode_zs, colorscale="ice", showscale=False, opacity=0.1
#         )
#         cathodes = draw_cathode_planes(
#             anode_xs, anode_ys, anode_zs, colorscale="burg", showscale=False, opacity=0.1
#         )
#         return center, anodes, cathodes


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


# def draw_light_detectors(data,event_id):
#     try:
#         charge = data["charge/events"][["id", "unix_ts"]][event_id]
#         num_light = data["light/events/data"].shape[0]
#         light = data["light/events", slice(0, num_light)][
#             ["id", "utime_ms"]
#         ]  # we have to try them all, events may not be time ordered
#     except:
#         print("No light information found, not plotting light detectors")
#         return []

#     match_light = match_light_to_charge_event(charge, light, event_id)

#     if match_light is None:
#         print(
#             f"No light event matches found for charge event {event_id}, not plotting light detectors"
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
#         charge = data["charge/events", self.file["id"]][["id", "unix_ts"]]
#         num_light = data["light/events/data"].shape[0]
#         light = data["light/events", slice(0, num_light)][
#             ["id", "utime_ms"]
#         ]  # we have to try them all, events may not be time ordered
#     except:
#         print("No light information found, not plotting light waveform")
#         return []

#     match_light = match_light_to_charge_event(charge, light, self.file["id"])

#     if match_light is None:
#         print(
#             f"No light event matches found for charge event {self.file['id']}, not plotting light waveform"
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
#     fig.update_yaxes(title_text='ADC counts')
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