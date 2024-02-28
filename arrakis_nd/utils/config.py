"""
Tools for parsing config files
"""
import yaml

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.utils import profiler


class ConfigParser:
    """
    """

    @profiler
    def __init__(
        self,
        config_file: str,
    ):
        """
        """
        self.config_file = config_file
        self.nested_config_files = [config_file]
        self.nested_configs = []
        self.data = {}
        self.parse_config()

    @profiler
    def parse_config(self):
        """
        """
        self.collect_nested_configs(self.config_file)
        self.nested_config_files.reverse()
        self.nested_configs.reverse()
        for config in self.nested_configs:
            self.data.update(config)

    @profiler
    def collect_nested_configs(
        self,
        config_file
    ):
        """_summary_

        Args:
            config_file (_type_): _description_
        """
        with open(config_file, 'r') as file:
            temp_data = yaml.safe_load(file)
            self.nested_configs.append(temp_data)
        if 'load_config' in temp_data.keys():
            self.nested_config_files.append(temp_data['load_config'])
            self.collect_nested_configs(temp_data['load_config'])

    @profiler
    def save_config(
        self,
        config_dictionary,
        output_file
    ):
        """_summary_

        Args:
            config_dictionary (_type_): _description_
            output_file (_type_): _description_
        """
        with open(output_file, 'w') as file:
            yaml.dump(config_dictionary, file)
