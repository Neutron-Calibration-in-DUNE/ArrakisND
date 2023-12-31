"""
ArrakisND main program
"""
import argparse
import torch

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.arrakis.arrakis import Arrakis


def run():
    """
    Arrkis main program.
    """
    parser = argparse.ArgumentParser(
        prog="Arrakis Module Runner",
        description="This program runs the Arrakis module " + "from a config file.",
        epilog="...",
    )
    parser.add_argument(
        "config_file",
        metavar="<str>.yml",
        type=str,
        help="config file specification for a BLIP module.",
    )
    parser.add_argument(
        "-name",
        dest="name",
        default="arrakis",
        help='name for this run (default "arrakis").',
    )
    parser.add_argument(
        "-n",
        dest="number_of_files",
        default=None,
        help='number of files to process (default None).',
    )
    args = parser.parse_args()
    # Setup config file.
    name = args.name
    number_of_files = args.number_of_files

    logger = Logger(name, output="both", file_mode="w")
    logger.info("configuring arrakis...")

    config = ConfigParser(args.config_file).data
    if "arrakis_nd" not in config.keys():
        logger.error("arrakis_nd section not in config!")
    arrakis_nd_config = config["arrakis_nd"]

    if number_of_files is not None:
        number_of_files = int(number_of_files)
        if isinstance(number_of_files, int):
            if number_of_files > 0:
                arrakis_nd_config["number_of_files"] = number_of_files

    system_info = logger.get_system_info()
    for key, value in system_info.items():
        logger.info(f"system_info - {key}: {value}")

    meta = {"config_file": args.config_file}
    if "verbose" in arrakis_nd_config:
        if not isinstance(arrakis_nd_config["verbose"], bool):
            logger.error(
                f'"arrakis_nd:verbose" must be of type bool, but got {type(arrakis_nd_config["verbose"])}!'
            )
        meta["verbose"] = arrakis_nd_config["verbose"]
    else:
        meta["verbose"] = False

    # Eventually we will want to check that the order of the arrakis_nds makes sense,
    # and that the data products are compatible and available for the different modes.

    # check for devices
    if "gpu" not in arrakis_nd_config.keys():
        logger.warn('"arrakis_nd:gpu" not specified in config!')
        gpu = None
    else:
        gpu = arrakis_nd_config["gpu"]
    if "gpu_device" not in arrakis_nd_config.keys():
        logger.warn('"arrakis_nd:gpu_device" not specified in config!')
        gpu_device = None
    else:
        gpu_device = arrakis_nd_config["gpu_device"]

    if torch.cuda.is_available():
        logger.info("CUDA is available with devices:")
        for ii in range(torch.cuda.device_count()):
            device_properties = torch.cuda.get_device_properties(ii)
            cuda_stats = f"name: {device_properties.name}, "
            cuda_stats += (
                f"compute: {device_properties.major}.{device_properties.minor}, "
            )
            cuda_stats += f"memory: {device_properties.total_memory}"
            logger.info(f" -- device: {ii} - " + cuda_stats)

    # set gpu settings
    if gpu:
        if torch.cuda.is_available():
            if gpu_device >= torch.cuda.device_count() or gpu_device < 0:
                logger.warn(
                    f"desired gpu_device '{gpu_device}' not available, using device '0'"
                )
                gpu_device = 0
            meta["device"] = torch.device(f"cuda:{gpu_device}")
            logger.info(
                f"CUDA is available, using device {gpu_device}"
                + f": {torch.cuda.get_device_name(gpu_device)}"
            )
        else:
            gpu is False
            logger.warn("CUDA not available! Using the cpu")
            meta["device"] = torch.device("cpu")
    else:
        logger.info("using cpu as device")
        meta["device"] = torch.device("cpu")
    meta["gpu"] = gpu

    arrakis = Arrakis(name, arrakis_nd_config, meta)
    arrakis.run_arrakis_nd()


if __name__ == "__main__":
    run()
