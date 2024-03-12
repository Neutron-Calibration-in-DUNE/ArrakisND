"""
This is an empty plugin that can be used as a template.
There are several comments throughout which provide detail
on the functions which must be overidden (user_name is replaced
with the $USER from the sesseion, but you should change it
to whatever you want).

name:   [*user_name*]
date:   [*user_date*]

"""
import h5py

from arrakis_nd.utils.utils import profiler
from arrakis_nd.plugins.plugin import Plugin


class EmptyPlugin(Plugin):
    """
    This custom plugin inherits from the base class
    Plugin, which is used in the main event loop of the
    Arrakis program.  Plugins are called by every worker
    that is processing events, and are called on every event.

    For each event, a pointer to the FLOW and ARRAKIS files
    are passed to the "process_event" function below.

    A dictionary of indices which point to the FLOW and ARRAKIS
    files is also passed as "event_indices" which contain the
    following keys:

    event_indices = {
        'interactions': indices of the "event" in the 'mc_truth/interactions/data' array
        'segments':     indices of the "event" in the 'mc_truth/segments/data' array
        'stack':        indices of the "event" in the 'mc_truth/stack/data' array
        'trajectories': indices of the "event" in the 'mc_truth/trajectories/data' array
        'charge':       indices of the "event" in the 'charge/calib_final_hits/data' array
        'light':        indices of the "event" in the 'light/sipm_hits/data' array
    }
    """
    def __init__(
        self,
        config: dict = {}
    ):
        """
        To keep things as simple as possible, we restrict plugins to
        have only a single input parameter "config", which is a dictionary
        that can contain any number of parameters that can be specified
        at run-time in the corresponding config file.

        The 'input_products' and 'output_produces' must be configured!
        The 'input_products' should be either 'None', 'str' or a list of 'str'.  The
        'str's tell Arrakis what plugin output_products that this plugin needs in
        order to work properly.  For example, the 'daughter_plugin' creates the
        'daughters' output_product, which contains the indices of the daughters
        for each traj_id in trajectories.  To specify that this plugin needs
        that product, we would write,

            self.input_products = 'daughters'   -  or  -  self.input_products = ['daughters']

        The 'output_products' should give the name of the intermediate product (if any)
        that this plugin produces and appends to the "event_products" dictionary.

        NB! Note that any changes one wishes to make to the arrakis_file happen
            in a subtle way.  For example, if you index any dataset from the arrakis_file,
            that will create a new numpy array and will NOT give you a reference
            to the original dataset.  Therefore, any changes you make must be sent back to
            the original arrakis_file in the end!  See the particle examples in the plugin
            folder.
        """
        super(EmptyPlugin, self).__init__(config)

        self.input_products = None
        self.output_products = None

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
        pass
