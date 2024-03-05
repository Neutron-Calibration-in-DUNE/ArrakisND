"""
"""
import h5py

from arrakis_nd.plugins.plugin import Plugin


class BarePlugin(Plugin):
    """
    An empty plugin to use as a starting point for a custom one.
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        """
        super(BarePlugin, self).__init__(config)

        self.input_products = None
        self.output_product = None

    def process_event(
        self,
        flow_file: h5py.File,
        arrakis_file: h5py.File,
        event_indices: dict,
        event_products: dict,
    ):
        """
        """
        pass
