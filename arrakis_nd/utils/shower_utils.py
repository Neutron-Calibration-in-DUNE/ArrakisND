"""
Utilities for fitting shower data
"""
import numpy as np
from scipy.interpolate import splprep, splev
import warnings
from scipy.integrate import quad, IntegrationWarning
from sklearn.decomposition import PCA

from arrakis_nd.utils.utils import profiler, integrand


@profiler
def fit_shower(
    shower_t0s,
    shower_xyz,
    ancestor_xyz_begin,
    ancestor_pxyz_begin
):
    """
    This function ...
    """
    try:
        """Normalize pxyz_begin"""
        ancestor_pxyz_begin_unit = ancestor_pxyz_begin / np.linalg.norm(ancestor_pxyz_begin)

        """Project charge onto the direction vector"""
        numerator = np.dot(shower_xyz - ancestor_xyz_begin, ancestor_pxyz_begin_unit)
        denominator = np.dot(ancestor_pxyz_begin_unit, ancestor_pxyz_begin_unit)
        t_values = numerator / denominator

        """Find the furthest point index"""
        furthest_point_index = np.argmax(np.abs(t_values))

        """Find the furthest point on the line"""
        furthest_point_on_line = ancestor_xyz_begin + t_values[furthest_point_index] * ancestor_pxyz_begin_unit

        """Calculate the shower length"""
        shower_length = np.sqrt(
            (ancestor_xyz_begin[0] - furthest_point_on_line[0])**2 +
            (ancestor_xyz_begin[1] - furthest_point_on_line[1])**2 +
            (ancestor_xyz_begin[2] - furthest_point_on_line[2])**2
        )

        """ Calculate the vector from p to each point"""
        difference_vectors = shower_xyz - ancestor_xyz_begin

        """ Calculate the cross product of x_unit with each vector"""
        cross_products = np.cross(difference_vectors, ancestor_pxyz_begin_unit)

        """ Calculate the Euclidean norm of each cross product vector"""
        """ This gives the area of the parallelogram formed by the vector and x_unit"""
        distances = np.linalg.norm(cross_products, axis=1)

        """ Divide by the magnitude of x to get the height of the parallelogram, which is the perpbeginicular distance"""
        perpbeginicular_distances = distances / np.linalg.norm(ancestor_pxyz_begin)

        """ Find the maximum perpbeginicular distance"""
        shower_width = np.max(perpbeginicular_distances)
        width_point = np.argmax(perpbeginicular_distances)
        width_point_on_line = ancestor_xyz_begin + t_values[width_point] * ancestor_pxyz_begin_unit
        width_distance = np.sqrt(
            (ancestor_xyz_begin[0] - width_point_on_line[0])**2 +
            (ancestor_xyz_begin[1] - width_point_on_line[1])**2 +
            (ancestor_xyz_begin[2] - width_point_on_line[2])**2
        )
        shower_angle = np.arctan(width_distance / shower_width)
        """Compute the shower angle"""

        """Sort charge with respect to t0s"""
        combined = sorted(zip(shower_t0s, shower_xyz), key=lambda x: x[0])
        sorted_t0, sorted_xyz = zip(*combined)
        sorted_t0 = np.array(sorted_t0)
        sorted_xyz = np.array(sorted_xyz)

        """Fill the values"""
        shower_data = {
            'shower_end': furthest_point_on_line,
            'shower_dir': [0, 0, 0],
            'len_gcm2': 0,
            'len_cm': shower_length,
            'width': shower_width,
            'angle': shower_angle,
        }
    except Exception:
        """Otherwise, these values are undefined"""
        shower_data = {
            'shower_end': [0, 0, 0],
            'shower_dir': [0, 0, 0],
            'len_gcm2': 0,
            'len_cm': 0,
            'width': 0,
            'angle': 0,
        }

    return shower_data
