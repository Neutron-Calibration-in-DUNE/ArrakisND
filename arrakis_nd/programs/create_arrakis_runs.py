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
    parser = argparse.ArgumentParser(
        prog='ArrakisND Config Generator',
        description='This program constructs a ...',
        epilog='...'
    )
    parser.add_argument(
        'config_file', metavar='<str>.yml', type=str,
        help='config file specification for a module.'
    )
    parser.add_argument(
        '-arrakis_config_location', dest='arrakis_config_location', default='/local_blip',
        help='location for the local blip directory.'
    )
    parser.add_argument(
        '-number_of_processes', dest='number_of_processes', default=1,
        help='number of processes for arrakis runs.'
    )

    logger = Logger('arrakis_generator', output="both", file_mode="w")

    args = parser.parse_args()
    config_file = args.config_file
    arrakis_config_location = args.arrakis_config_location
    number_of_processes = int(args.number_of_processes)

    config_parser = ConfigParser(config_file)
    config = config_parser.data

    arrakis_dict = config["arrakis_nd"]
    simulation_folder = arrakis_dict["simulation_folder"]
    simulation_files = arrakis_dict["simulation_files"]
    if isinstance(arrakis_dict["simulation_files"], list):
        simulation_files = [
            input_file for input_file in arrakis_dict["simulation_files"]
            if input_file not in arrakis_dict["skip_files"]
        ]
    elif isinstance(arrakis_dict["simulation_files"], str):
        if arrakis_dict["simulation_files"] == "all":
            logger.info(
                f"searching {simulation_folder} recursively for all .h5 files."
            )
            simulation_files = [
                os.path.basename(input_file) for input_file in glob.glob(
                    f"{simulation_folder}*.h5", recursive=True
                )
                if input_file not in arrakis_dict["skip_files"]
            ]
        else:
            try:
                logger.info(
                    f'searching {simulation_folder} recursively for all {arrakis_dict["simulation_files"]} files.'
                )
                simulation_files = [
                    os.path.basename(input_file) for input_file in glob.glob(
                        f'{simulation_folder}{arrakis_dict["simulation_files"]}',
                        recursive=True,
                    )
                    if input_file not in arrakis_dict["skip_files"]
                ]
            except Exception as exception:
                logger.error(
                    f'specified "simulation_files" parameter: {arrakis_dict["simulation_files"]} incompatible!'
                    + f" exception: {exception}"
                )
    else:
        logger.error(
            f'specified "simulation_files" parameter: {arrakis_dict["simulation_files"]} incompatible!'
        )
    # split up simulation_files
    file_chunk_size = int(np.ceil(len(simulation_files) / number_of_processes))
    file_chunks = [simulation_files[i:i+file_chunk_size] for i in range(0, len(simulation_files), file_chunk_size)]

    arrakis_folders = []

    for ii in range(len(file_chunks)):
        random_config = config.copy()
        random_config['arrakis_nd']['simulation_files'] = file_chunks[ii]
        if not os.path.isdir(f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}'):
            os.makedirs(f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}')
        arrakis_folders.append([f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}'])
        config_parser.save_config(
            random_config,
            f'{os.path.abspath(arrakis_config_location)}/arrakis_config_iteration_{ii}/arrakis_config.yaml'
        )
    with open(f'{os.path.abspath(arrakis_config_location)}/arrakis_data.csv', "w") as file:
        writer = csv.writer(file, delimiter=",")
        writer.writerows(arrakis_folders)


if __name__ == "__main__":
    run()
