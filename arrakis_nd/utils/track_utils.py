"""
Utilities for fitting track data
"""
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning

from arrakis_nd.utils.utils import profiler, integrand


@profiler
def fit_track(
    track_t0s,
    track_xyz
):
    """
    This function ...
    """
    try:
        """Sort charge with respect to t0s"""
        combined = sorted(zip(track_t0s, track_xyz), key=lambda x: x[0])
        sorted_t0, sorted_xyz = zip(*combined)
        sorted_t0 = np.array(sorted_t0)
        sorted_xyz = np.array(sorted_xyz)

        """Try to fit a spline curve to the xyz data"""
        tck, u = splprep(
            [sorted_xyz[:, 0], sorted_xyz[:, 1], sorted_xyz[:, 2]],
            s=0
        )

        """Try to integrate this curve"""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", IntegrationWarning)
            try:
                curve_length, _ = quad(integrand, 0, 1, args=(tck,))
            except Exception:
                curve_length = 0

        """Get the derivatives at the beginning and ending points"""
        dxdt_start, dydt_start, dzdt_start = splev(0, tck, der=1)
        dxdt_end, dydt_end, dzdt_end = splev(1, tck, der=1)
        dx_start = np.array([dxdt_start, dydt_start, dzdt_start])
        dx_end = np.array([dxdt_end, dydt_end, dzdt_end])
        dx_start_magnitude = np.linalg.norm(dx_start)
        dx_end_magnitude = np.linalg.norm(dx_end)

        if dx_start_magnitude > 0:
            track_dir = [dxdt_start, dydt_start, dzdt_start] / dx_start_magnitude
        else:
            track_dir = [dxdt_start, dydt_start, dzdt_start]

        if dx_end_magnitude > 0:
            track_enddir = [dxdt_end, dydt_end, dzdt_end] / dx_end_magnitude
        else:
            track_enddir = [dxdt_end, dydt_end, dzdt_end]

        """Determine the integrated columnar density"""
        lar_rho = 1.40
        curve_length_rho = lar_rho * curve_length

        """Fill the values"""
        track_data = {
            'track_dir': track_dir,
            'track_enddir': track_enddir,
            'track_len_gcm2': curve_length_rho,
            'track_len_cm': curve_length
        }
    except Exception:
        """Otherwise, these values are undefined"""
        track_data = {
            'track_dir': [0, 0, 0],
            'track_enddir': [0, 0, 0],
            'track_len_gcm2': 0,
            'track_len_cm': 0
        }

    return track_data
