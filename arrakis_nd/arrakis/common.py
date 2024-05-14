"""
"""
import numpy as np

"""
Second, we set up arrays for CAF-like reco objects which consist
of the following types:
    (1) tracklettes - pieces of contiguous track like objects.
        (a) start (x,y,z) - geant4 track start point
        (b) end (x,y,z) - geant4 track end point
        (c) start_hit (x,y,z) - corresponding hit start point
        (d) end_hit (x,y,z) - corresponding hit end point
        (e) dir (u_x,u_y,u_z) - geant4 track momentum start dir
        (f) enddir (u_x,u_y,u_z) - geant4 track momentum end dir
        (g) dir_hit (u_x,u_y,u_z) - estimate of track direction taken from start_hit point
        (h) enddir_hit (u_x,u_y,u_z) - estimate of track direction takend from end_hit point
        (i) Evis - visible energy in voxels corresponding to this track
        (j) qual - reco specific quality metric
        (k) len_gcm2 - track length in g/cm2
        (l) len_cm - track length in centimeter
        (m) E - track energy estimate in MeV
        (n) truth - associated true particle in flow file (if relevant) (i.e. track_id)
        (o) truthOverlap - fractional overlap between this track and true particle
    (2) tracks - collections of tracklettes which define a complete track.
    (3) fragments - pieces of contiguous shower like objects.
        (a) start (x,y,z) - geant4 shower start point
        (b) start_hit (x,y,z) - corresponding hit start point
        (c) dir (u_x,u_y,u_z) - geant4 shower direction
        (d) dir_hit (u_x,u_y,u_z) - estimate of shower direction taken from start_hit point
        (e) Evis - visible energy in voxels corresponding to this shower
        (f) qual - reco specific quality metric
        (g) len_gcm2 - track length in g/cm2
        (h) len_cm - track length in centimeter
        (i) E - track energy estimate in MeV
        (j) truth - associated true particle in flow file (if relevant) (i.e. track_id)
        (k) truthOverlap - fractional overlap between this shower and true particle
    (4) showers - collections of fragments which define a complete shower.
    (5) blips - collections of points associated to blip-like activity.
        (a) start (x,y,z) - geant4 start position of this blip object
        (b) start_hit (x,y,z) - corresponding hit start point
        (c) Evis - visible energy in voxels corresponding to this blip
        (d) E - reconstructed energy (GeV)
        (e) bliphyp - hypothesis for this blip's identity
        (f) truth - associated true particle in flow file (if relevant) (i.e track_id)
        (g) truthOverlap - fractional overlap between this blip and true particle
    (6) particles - particles within an event that are associated to tracks/showers/blips.
        (a) primary - is this reco particle a primary one (i.e. eminates directly from vertex)?
        (b) pdg - pdg code inferred for this particle
        (c) tgtA - atomic number of nucleus this particle was reconstructed in
        (d) score - PID score for this particle
        (e) Evis - visible energy in voxels corresponding to this particle
        (f) E - reconstructed energy (GeV)
        (g) E_method - method used to determine energy for the particle
        (h) p (p_x,p_y,p_z) - reconstructed momentum for this particle
        (i) start (x,y,z) - reconstructed start point of this particle
        (j) end (x,y,z) - reconstructed end point of this particle
        (k) start_hit (x,y,z) - reconstructed start hit of this particle
        (l) end_hit (x,y,z) - reconstructed end hit of this particle
        (m) contained - contained in LAr TPC?
        (n) truth - associated true particle in flow file (if relevant) (i.e. track_id)
        (o) truthOverlap - fractional overlap between this reco particle and true particle
    (7) interactions - collections of tracks/showers/blips which define interactions of interest.
        (a) vtx (x,y,z) - reconstructed vertex location
        (b) dir (u_x,u_y,u_z) - hypothesis for this interaction's parent particle direction
            (b.i) lngtrk - direction using longest track
            (b.ii) heshw - direction using highest energy shower
        (c) nuhyp - hypothesis for this interaction's neutrino identity
            (c.i) isnubar
            (c.ii) nue
            (c.iii) numu
            (c.iv) nutau
            (c.v) nc
            (c.vi) protons0
            (c.vii) protons1
            (c.viii) protons2
            (c.ix) protonsN
            (c.x) chgpi0
            (c.xi) chgpi1
            (c.xii) chgpi2
            (c.xiii) chgpiN
            (c.xiv) pizero0
            (c.xv) pizero1
            (c.xvi) pizero2
            (c.xvii) pizeroN
            (c.xviii) neutron0
            (c.xix) neutron1
            (c.xx) neutron2
            (c.xxi) neutronN
        (d) Enu - hypothesis for this interaction's neutrino energy
        (e) part - collections of reconstructed particles
        (f) truth - indicies of true_interactions in flow file (if relevant)
        (g) truthOverlap - fractional overlap between this reco interaction and each true interaction
"""

tracklette_data_type = np.dtype([
    ('event_id', 'i4'),
    ('tpc_id', 'i4'),
    ('tracklette_id', 'i4'),
    ('vertex_id', 'i8'),
    ('tracklette_type', 'i4'),
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

track_data_type = np.dtype([
    ('event_id', 'i4'),
    ('track_id', 'i4'),
    ('vertex_id', 'i8'),
    ('track_type', 'i4'),
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

fragment_data_type = np.dtype([
    ('event_id', 'i4'),
    ('tpc_id', 'i4'),
    ('fragment_id', 'i4'),
    ('vertex_id', 'i8'),
    ('ancestor_id', 'i4'),
    ('fragment_type', 'i4'),
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

shower_data_type = np.dtype([
    ('event_id', 'i4'),
    ('shower_id', 'i4'),
    ('vertex_id', 'i8'),
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
    ('width', 'f4'),
    ('angle', 'f4'),
    ('E', 'f4'),
    ('truth', 'i4', (1, 20)),
    ('truthOverlap', 'f4', (1, 20)),
])

blip_data_type = np.dtype([
    ('event_id', 'i4'),
    ('tpc_id', 'i4'),
    ('blip_id', 'i4'),
    ('vertex_id', 'i8'),
    ('ancestor_id', 'i4'),
    ('blip_type', 'i4'),
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

particle_data_type = np.dtype([
    ('event_id', 'i4'),
    ('particle_id', 'i4'),
    ('vertex_id', 'i8'),
    ('primary', 'i4'),
    ('pdg', 'i4'),
    ('tgtA', 'i4'),
    ('score', 'f4'),
    ('Evis', 'f4'),
    ('E', 'f4'),
    ('E_method', 'i4'),
    ('p', 'f4', (1, 3)),
    ('dir_hit', 'f4', (1, 3)),
    ('start', 'f4', (1, 3)),
    ('end', 'f4', (1, 3)),
    ('start_hit', 'f4', (1, 3)),
    ('end_hit', 'f4', (1, 3)),
    ('contained', 'i4'),
    ('truth', 'i4', (1, 20)),
    ('truthOverlap', 'f4', (1, 20)),
])

neutrino_data_type = np.dtype([
    ('event_id', 'i4'),
    ('neutrino_id', 'i4'),
    ('vtx', 'f4', (1, 3)),
    ('dir_lngtrk', 'f4', (1, 3)),
    ('dir_heshw', 'f4', (1, 3)),
    ('nuhyp', 'f4', (1, 20)),
    ('Enu', 'f4'),
    ('E_method', 'i4'),
    ('part', 'i4', (1, 20)),
    ('truth', 'i4', (1, 20)),
    ('truthOverlap', 'f4', (1, 20)),
])

interaction_data_type = np.dtype([
    ('event_id', 'i4'),
    ('interaction_id', 'i4'),
    ('vtx', 'f4', (1, 3)),
    ('dir_lngtrk', 'f4', (1, 3)),
    ('dir_heshw', 'f4', (1, 3)),
    ('nuhyp', 'f4', (1, 20)),
    ('Enu', 'f4'),
    ('E_method', 'i4'),
    ('part', 'i4', (1, 20)),
    ('truth', 'i4', (1, 20)),
    ('truthOverlap', 'f4', (1, 20)),
])

neutrino_event_data_type = np.dtype([
    ('event_id', 'i4'),
    ('neutrino_type', 'i4'),
    ('fiducialized', 'i4')
])

interaction_event_data_type = np.dtype([
    ('event_id', 'i4'),
    ('interaction_type', 'i4'),
    ('fiducialized', 'i4')
])

nar_inelastic_data_type = np.dtype([
    ('event_id', 'i4'),
    ('vertex_id', 'i8'),
    ('proton_id', 'i4'),
    ('proton_xyz_start', 'f4', (1, 3)),
    ('nu_vertex', 'f4', (1, 3)),
    ('proton_total_energy', 'f4'),
    ('proton_vis_energy', 'f4'),
    ('proton_length', 'f4'),
    ('proton_pxyz_start', 'f4', (1, 3)),
    ('proton_pxyz_end', 'f4', (1, 3)),
    ('parent_total_energy', 'f4'),
    ('parent_length', 'f4'),
    ('parent_pxyz_start', 'f4', (1, 3)),
    ('parent_pxyz_end', 'f4', (1, 3)),
    ('nu_proton_dt', 'f4'),
    ('nu_proton_distance', 'f4'),
    ('parent_pdg', 'i4'),
    ('grandparent_pdg', 'i4'),
    ('primary_pdg', 'i4'),
    ('primary_length', 'f4'),
    ('neutrino_vertex_fiducialized', 'i4'),
    ('proton_xyz_start_fiducialized', 'i4'),
    ('proton_xyz_end_fiducialized', 'i4'),
    ('truth_segment_overlap', 'f4'),
    ('best_completeness_cluster', 'i4'),
    ('best_completeness', 'f4'),
    ('best_purity', 'f4')
])

undefined_data_type = np.dtype([
    ('event_id', 'i4'),
    ('vertex_id', 'i8'),
    ('traj_id', 'i8'),
    ('parent_id', 'i8'),
    ('pdg_id', 'i8'),
    ('parent_pdg_id', 'i8'),
    ('ancestor_id', 'i8'),
    ('ancestor_pdg_id', 'i8'),
    ('ancestor_level', 'i4'),
    ('start_process', 'i4'),
    ('start_subprocess', 'i4'),
    ('end_process', 'i4'),
    ('end_subprocess', 'i4'),
    ('parent_start_process', 'i4'),
    ('parent_start_subprocess', 'i4'),
    ('parent_end_process', 'i4'),
    ('parent_end_subprocess', 'i4')
])

output_error_data_type = np.dtype([
    ('event_id', 'i4'),
    ('plugin_name', 'S50'),
    ('plugin_error', 'S500'),
    ('plugin_line_number', 'S50'),
    ('plugin_file_name', 'S50'),
])