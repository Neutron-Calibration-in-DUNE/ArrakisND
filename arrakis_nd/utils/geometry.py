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
            self.active_tpc_widths_2x2 = []

            self.parse_config()

        
        def parse_config(self):
            self.active_tpc_widths = self.active_tpc_widths_2x2
        
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