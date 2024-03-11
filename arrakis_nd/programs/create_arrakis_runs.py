"""
Script for generating hyperparameter config files from a base config
"""
import os
import numpy as np
import argparse
import csv
import glob

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser


def run():
    """
    This program splits up an Arrakis config over several config
    files by distributing the flow_files into separate jobs.  This
    allows Arrakis to run on multiple files at once, which greatly
    speeds up runtime.
    """
    parser = argparse.ArgumentParser(
        prog='ArrakisND Config Generator',
        description='This program splits up an Arrakis config over several config \
            files by distributing the flow_files into separate jobs.  This \
            allows Arrakis to run on multiple files at once, which greatly \
            speeds up runtime.',
        epilog='...'
    )
    parser.add_argument(
        'config_file',
        metavar='<str>.yml',
        type=str,
        help='config file specification for a module.'
    )
    parser.add_argument(
        '-arrakis_config_location',
        dest='arrakis_config_location',
        default='/local_arrakis',
        help='location for the local arrakis directory.'
    )
    parser.add_argument(
        '-number_of_processes',
        dest='number_of_processes',
        default=1,
        help='number of processes for arrakis runs.'
    )

    """Set up logger"""
    logger = Logger('arrakis_generator')

    """Parse command line arguments"""
    args = parser.parse_args()
    config_file = args.config_file
    arrakis_config_location = args.arrakis_config_location
    number_of_processes = int(args.number_of_processes)

    """Pull in config file"""
    config_parser = ConfigParser(config_file)
    config = config_parser.data

    arrakis_dict = config["arrakis_nd"]
    flow_folder = arrakis_dict["flow_folder"]
    flow_files = arrakis_dict["flow_files"]

    if isinstance(arrakis_dict["flow_files"], list):
        """
        If the flow_files parameter is a list, look through
        the list and make sure each specified file actually exists
        in the flow_folder.
        """
        flow_files = [
            input_file for input_file in arrakis_dict["flow_files"]
            if input_file not in arrakis_dict["skip_files"]
        ]
        for flow_file in flow_files:
            if not os.path.isfile(flow_file):
                logger.error(
                    f"specified file {flow_file} does not exist in directory {flow_folder}!"
                )
    elif isinstance(arrakis_dict["flow_files"], str):
        """
        If the flow_files parameter is a string, check if its
        the phrase 'all', and if so, recursively grab all h5
        files in the flow_folder.

        Otherwise, assume that the flow_files parameter is a
        file extension, and search recursively for all files
        with that extension.
        """
        if arrakis_dict["flow_files"] == "all":
            logger.info(
                f"searching {flow_folder} recursively for all .h5 FLOW files."
            )
            flow_files = [
                os.path.basename(input_file) for input_file in glob.glob(
                    f"{flow_folder}*.h5", recursive=True
                )
                if 'FLOW' in input_file and input_file not in arrakis_dict["skip_files"]
            ]
        else:
            try:
                logger.info(
                    f'searching {flow_folder} recursively for all {arrakis_dict["flow_files"]} files.'
                )
                flow_files = [
                    os.path.basename(input_file) for input_file in glob.glob(
                        f'{flow_folder}/{arrakis_dict["flow_files"]}',
                        recursive=True,
                    )
                    if input_file not in arrakis_dict["skip_files"]
                ]
            except Exception as exception:
                logger.error(
                    f'specified "flow_files" parameter: {arrakis_dict["flow_files"]} incompatible!'
                    + f" exception: {exception}"
                )
    else:
        logger.error(
            f'specified "flow_files" parameter: {arrakis_dict["flow_files"]} incompatible!'
        )

    """Split up flow files among the jobs"""
    file_chunk_size = int(np.ceil(len(flow_files) / number_of_processes))
    file_chunks = [flow_files[i:i+file_chunk_size] for i in range(0, len(flow_files), file_chunk_size)]

    arrakis_folders = []

    """Generate the corresponding config files for each job"""
    for ii in range(len(file_chunks)):
        random_config = config.copy()
        random_config['arrakis_nd']['flow_files'] = file_chunks[ii]
        if not os.path.isdir(f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}'):
            os.makedirs(f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}')
        arrakis_folders.append([f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}'])
        config_parser.save_config(
            random_config,
            f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}/arrakis_config.yaml'
        )

    """Save config file information so that the job distribution script knows where to look"""
    with open(f'{os.path.abspath(arrakis_config_location)}/arrakis_data.csv', "w") as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerows(arrakis_folders)


if __name__ == "__main__":
    run()
