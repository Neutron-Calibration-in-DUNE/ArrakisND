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

    arrakis_dict = config['arrakis_nd']

    """Check for parameters"""
    if 'flow_folder' not in arrakis_dict.keys():
        logger.error('flow_folder not specified in config!')
    if 'flow_files' not in arrakis_dict.keys():
        logger.error('flow_files not specified in config!')

    flow_folder = arrakis_dict['flow_folder']
    flow_files = arrakis_dict["flow_files"]

    """Check for arrakis folder"""
    if 'arrakis_folder' not in arrakis_dict.keys():
        logger.warn('arrakis_folder not specified in config! setting to "/local_scratch"')
        arrakis_dict['arrakis_folder'] = '/local_scratch'
    arrakis_folder = arrakis_dict['arrakis_folder'].replace('flow', 'arrakis').replace('FLOW', 'ARRAKIS')

    """Check that flow folder exists"""
    if not os.path.isdir(flow_folder):
        logger.error(f'specified flow_folder {flow_folder} does not exist!')

    """Check that flow folder has a '/' at the end"""
    if flow_folder[-1] != '/':
        flow_folder += '/'

    """Check that arrakis folder has a '/' at the end"""
    if arrakis_folder[-1] != '/':
        arrakis_folder += '/'

    """Check that arrakis folder exists"""
    if not os.path.isdir(arrakis_folder):
        os.makedirs(arrakis_folder)

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
                        f'{flow_folder}/*.{arrakis_dict["flow_files"]}',
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
    if len(flow_files) == 0:
        logger.error(f"found no flow files in {flow_folder}!")
    logger.info(f"found {len(flow_files)} flow files")

    """Split up flow files among the jobs"""
    file_chunk_size = len(flow_files) // number_of_processes
    larger_file_chunk_size = len(flow_files) % number_of_processes
    logger.info(f"file_chunk_size is {file_chunk_size}")
    logger.info(f"remainder size is {larger_file_chunk_size}")

    file_chunks = []
    start = 0
    for ii in range(number_of_processes):
        end = start + file_chunk_size + (1 if ii < larger_file_chunk_size else 0)
        file_chunks.append(flow_files[start:end])
        start = end
    logger.info(f"split up {len(flow_files)} into {number_of_processes}")

    """Generate the corresponding config files for each job"""
    for ii in range(len(file_chunks)):
        config_iteration = config.copy()
        config_iteration['arrakis_nd']['flow_files'] = file_chunks[ii]
        config_parser.save_config(
            config_iteration,
            f'{os.path.abspath(arrakis_config_location)}/arrakis_config_{ii}.yaml'
        )
    logger.info(f"saving split configs into {arrakis_config_location}")
    logger.info("create_arrakis_runs job complete. closing out")


if __name__ == "__main__":
    run()
