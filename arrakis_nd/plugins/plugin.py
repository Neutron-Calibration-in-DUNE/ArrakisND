"""
"""
import h5py


class Plugin:
    """
    A base class for plugins to Arrakis.

    Plugins are classes which process flow data and produce
    either intermediate data products or arrakis output.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        All plugins must run super(CustomPlugin, self).__init__(config)
        and specify any required data_products and any output data_products
        so that Arrakis knows they are called in the correct sequence.

        Args:
            config (dict, optional): _description_. Defaults to {}.
        """
        self.config = config

        """These should be set by the creator of the custom plugin"""
        self.input_products = None
        self.output_product = None

    def process_event(
        self,
        event: int,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        Event processing function for a plugin. The process_event
        function requires both a flow file and arrakis file, as
        well as the dictionary of event_indices corresponding to
        both.  The event_products is a shared dictionary which
        can be updated with intermmediate data used by other
        plugins.


        This function must be overriden by the inherited class.

        Args:
            flow_file (h5py.File): _description_
            arrakis_file (h5py.File): _description_
            event_indices (dict): _description_
            event_products (dict): _description_
        """
        raise NotImplementedError("Plugin classes must implement this function!")
