"""
ArrakisND main program
"""
import argparse
import torch

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.arrakis.arrakis import Arrakis


def run():
    """_summary_
    """
    """
    Arrakis main program.
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

    config = ConfigParser(args.config_file).data
    if number_of_files is not None:
        number_of_files = int(number_of_files)
        if isinstance(number_of_files, int):
            if number_of_files > 0:
                config["arrakis_nd"]["number_of_files"] = number_of_files
    meta = {
        "name": name,
        "config_file": args.config_file
    }

    arrakis = Arrakis(config, meta)
    arrakis.run_arrakis_nd()


if __name__ == "__main__":
    """_summary_
    """
    run()
