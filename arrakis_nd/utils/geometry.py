    """
    """
    import numpy as np
    
    
    class Geometry:
        """
        This class contains several useful functions related
        to the detector geometry of the various ND TPCs. 
        """
        def __init__(
            self,
            config: dict = {}
        ):
            self.config = config

            """Configuration for 2x2 geometry"""
            self.active_tpc_widths_2x2 = [30.6, 130., 64.]
            self.tpcs_relative_to_module_2x2=[
                [-15.7, 0., 0.],
                [15.7,  0., 0.]
            ]
            self.modules_relative_to_2x2=[
                [-33.5,0., -33.5],
                [33.5,0.,  -33.5],
                [-33.5,0., 33.5],
                [33.5,0.,  33.5]
            ]
            self.detector_center_2x2=[0., -268, 1300]
            
            """Set up detector geometry"""
            self.parse_config()

        def parse_config(self):
            """
            Here we determine what the various detector variables should be,
            and then construct other variables relative to them.
            """
            self.active_tpc_widths = self.active_tpc_widths_2x2
            self.tpcs_relative_to_module = self.tpcs_relative_to_module_2x2
            self.modules_relative_to = self.modules_relative_to_2x2
            self.detector_center = self.detector_center_2x2
        
        def determine_tpc_bounds(i):
            active_tpc_widths=[30.6, 130., 64.] # cm
            tpcs_relative_to_module=[[-15.7,0.,0.],[15.7,0.,0.]]
            modules_relative_to_2x2=[[-33.5,0.,-33.5],
                                    [33.5,0.,-33.5],
                                    [-33.5,0.,33.5],
                                    [33.5,0.,33.5]]
            detector_center=[0.,-268,1300]
            tpc_bounds=np.array([-active_tpc_widths[i]/2., active_tpc_widths[i]/2.])
            tpc_bounds_relative_to_2x2=[]
            for tpc in tpcs_relative_to_module:
                tpc_bound_relative_to_module = tpc_bounds + tpc[i]
                for module in modules_relative_to_2x2:
                    bound = tpc_bound_relative_to_module + module[i]
                    tpc_bounds_relative_to_2x2.append(bound)
            bounds_relative_to_NDhall = np.array(tpc_bounds_relative_to_2x2) + detector_center[i]
            return np.unique(bounds_relative_to_NDhall, axis=0)