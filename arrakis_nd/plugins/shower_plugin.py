"""
"""
import h5py
import numpy as np

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin
from arrakis_nd.arrakis.common import shower_data_type
from arrakis_nd.dataset.common import (
    ProcessType, SubProcessType,
    Shower
)


class ShowerPlugin(Plugin):
    """_summary_

    Args:
        Plugin (_type_): _description_
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(ShowerPlugin, self).__init__(config)

        self.input_products = [
            'daughters',
            'track_id_hit_map',
            'track_id_hit_segment_map',
            'track_id_hit_t0_map',
            'parent_pdg_id',
            'ancestor_traj_id_map',
            'ancestor_pdg_id_map',
        ]
        self.output_products = ['shower']

        if "shower_size_threshold" not in self.config:
            self.config["shower_size_threshold"] = 10
        self.shower_size_threshold = self.config["shower_size_threshold"]

    @profiler
    def process_event(
        self,
        event: int,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        """
        trajectories = flow_file['mc_truth/trajectories/data'][event_indices['trajectories']]

        trajectories_traj_ids = trajectories['traj_id']
        trajectories_vertex_ids = trajectories['vertex_id']
        trajectories_pdg_ids = trajectories['pdg_id']
        trajectories_xyz_start = trajectories['xyz_start']
        trajectories_xyz_end = trajectories['xyz_end']
        trajectories_pxyz_start = trajectories['pxyz_start']
        trajectories_pxyz_end = trajectories['pxyz_end']
        trajectories_E = trajectories['E_start']
        trajectories_start_subprocess = trajectories['start_subprocess']

        """Collect unique fragments"""
        arrakis_fragments = event_products['fragment']
        arrakis_blips = event_products['blip']
        ancestor_traj_index_map = event_products['ancestor_traj_index_map']
        ancestor_traj_id_map = event_products['ancestor_traj_id_map']
        ancestor_pdg_id_map = event_products['ancestor_pdg_id_map']

        fragment_ids = np.array([fragment['fragment_id'][0] for fragment in arrakis_fragments])
        fragment_vertex_ids = np.array([fragment['vertex_id'][0] for fragment in arrakis_fragments])
        blip_ids = np.array([blip['blip_id'][0] for blip in arrakis_blips])
        blip_vertex_ids = np.array([blip['vertex_id'][0] for blip in arrakis_blips])

        fragment_ancestors = {}
        ancestor_vertex = {}
        blip_ancestors = {}

        for ii, (fragment_id, vertex_id) in enumerate(zip(
            fragment_ids, fragment_vertex_ids
        )):
            ancestor_traj_index = ancestor_traj_index_map[(fragment_id, vertex_id)]
            ancestor_id = ancestor_traj_id_map[(fragment_id, vertex_id)]
            ancestor_vertex[ancestor_id] = vertex_id
            if ancestor_id not in fragment_ancestors.keys():
                fragment_ancestors[ancestor_id] = [ii]
                blip_ancestors[ancestor_id] = []
            else:
                fragment_ancestors[ancestor_id].append(ii)

        for ii, (blip_id, vertex_id) in enumerate(zip(
            blip_ids, blip_vertex_ids
        )):
            ancestor_traj_index = ancestor_traj_index_map[(blip_id, vertex_id)]
            ancestor_id = ancestor_traj_id_map[(blip_id, vertex_id)]
            if ancestor_id in blip_ancestors.keys():
                blip_ancestors[ancestor_id].append(ii)

        """Now that we've collected unique ancestors, lets make a shower for each"""
        for ii, (ancestor_id, frag_ids) in enumerate(fragment_ancestors.items()):
            ancestor_traj_index = ancestor_traj_index_map[(ancestor_id, ancestor_vertex[ancestor_id])]
            ancestor_pdg = trajectories_pdg_ids[ancestor_traj_index]
            ancestor_xyz_start = trajectories_xyz_start[ancestor_traj_index]
            ancestor_pxyz_start = trajectories_pxyz_start[ancestor_traj_index]

            """First we check to see whether a conversion exists, if not then no shower"""
            num_conversions = 0
            for frag_id in frag_ids:
                traj_index = np.where(
                    (trajectories_traj_ids == arrakis_fragments[frag_id]['fragment_id']) &
                    (trajectories_vertex_ids == arrakis_fragments[frag_id]['vertex_id'])
                )
                if trajectories_start_subprocess[traj_index] == SubProcessType.GammaConversion.value:
                    num_conversions += 1

            """Determine shower type"""
            if (
                (abs(ancestor_pdg) == 11) |
                (abs(ancestor_pdg) == 22)
            ):
                shower_type = Shower.Electromagnetic.value
            elif ancestor_pdg == 111:
                shower_type = Shower.Pion0Decay.value
            elif ancestor_pdg == 211:
                shower_type = Shower.PiPlusDecay.value
            elif ancestor_pdg == -211:
                shower_type = Shower.PiMinusDecay.value
            elif ancestor_pdg == 311:
                shower_type = Shower.Kaon0Decay.value
            elif ancestor_pdg == 130:
                shower_type = Shower.KaonLongDecay.value
            elif ancestor_pdg == 310:
                shower_type = Shower.KaonShortDecay.value
            elif ancestor_pdg == 321:
                shower_type = Shower.KaonPlusDecay.value
            elif ancestor_pdg == -321:
                shower_type = Shower.KaonMinusDecay.value
            elif ancestor_pdg == 421:
                shower_type = Shower.D0Decay.value
            elif ancestor_pdg == 411:
                shower_type = Shower.DPlusDecay.value
            elif ancestor_pdg == -411:
                shower_type = Shower.DMinusDecay.value
            elif ancestor_pdg == 3122:
                shower_type = Shower.LambdaDecay.value
            elif ancestor_pdg == 3212:
                shower_type = Shower.Sigma0Decay.value
            elif ancestor_pdg == 3222:
                shower_type = Shower.SigmaPlusDecay.value
            elif ancestor_pdg == 3112:
                shower_type = Shower.SigmaMinusDecay.value
            else:
                shower_type = Shower.Undefined.value

            """Determine start and end of shower"""
            """
            Create the shower object
                shower_data_type = np.dtype([
                    ('event_id', 'i4'),
                    ('shower_id', 'i4'),
                    ('vertex_id', 'i8'),
                    ('blip_ids', 'i4', (1, 400)),
                    ('shower_type', 'i4'),
                    ('start', 'f4', (1, 3)),
                    ('end', 'f4', (1, 3)),
                    ('start_hit', 'f4', (1, 3)),
                    ('end_hit', 'f4', (1, 3)),
                    ('dir', 'f4', (1, 3)),
                    ('enddir', 'f4', (1, 3)),
                    ('dir_hit', 'f4', (1, 3)),
                    ('enddir_hit', 'f4', (1, 3)),
                    ('Evis', 'f4'),
                    ('qual', 'f4'),
                    ('len_gcm2', 'f4'),
                    ('len_cm', 'f4'),
                    ('E', 'f4'),
                    ('truth', 'i4', (1, 20)),
                    ('truthOverlap', 'f4', (1, 20)),
                ])
            """

            shower_data = np.array([(
                event,
                ancestor_id,
                vertex_id,
                shower_type,
                ancestor_xyz_start,
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                ancestor_pxyz_start,
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
                0,
                1,
                0,
                0,
                0,
                [[ancestor_traj_index]+[0]*19],
                [[1]+[0]*19])],
                dtype=shower_data_type
            )

            """Add the new shower to the CAF objects"""
            event_products['shower'].append(shower_data)
