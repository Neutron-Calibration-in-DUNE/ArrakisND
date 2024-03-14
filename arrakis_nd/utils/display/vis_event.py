import plotly.graph_objects as go
import dash_core_components as dcc
import dash_html_components as html

def create_3d_scatter():
    # Create a 3D scatter plot
    fig = go.Figure(data=[go.Scatter3d(
        x=[1, 2, 3],
        y=[4, 5, 6],
        z=[7, 8, 9],
        mode='markers',
        marker=dict(
            size=12,
            color=[1, 2, 3],  # set color to an array/list of desired values
            colorscale='Viridis',  # choose a colorscale
            opacity=0.8
        )
    )])

    # Set the layout of the plot
    fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))

    # Return a Div containing the plot
    return html.Div([
        html.P("3D DISPLAY OF THE SCATTER PLOT:"),
        dcc.Graph(
            id='3d-scatter',
            figure=fig
        )
    ])