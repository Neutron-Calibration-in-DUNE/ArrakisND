# interactive_dashboard.py

import panel as pn
import holoviews as hv
import numpy as np

# Ensure compatibility for running in notebooks
hv.extension('bokeh')
pn.extension()

# Define your interactive function
def interactive_scatter(n_points=100):
    x = np.random.randn(n_points)
    y = np.random.randn(n_points)
    return hv.Scatter((x, y)).opts(size=5, color='blue')

# Create a slider widget
slider = pn.widgets.IntSlider(name='Number of points', start=10, end=1000, step=10, value=100)

# Bind the function and slider
scatter_plot = pn.bind(interactive_scatter, n_points=slider)

# Define the layout
layout = pn.Column(slider, scatter_plot)