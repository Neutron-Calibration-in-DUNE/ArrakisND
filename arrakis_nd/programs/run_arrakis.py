"""
ArrakisND main program
"""
import argparse
from mpi4py import MPI

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.config import ConfigParser
from arrakis_nd.arrakis.arrakis import Arrakis


def run():
    """_summary_
    """

    """
    We do a preliminary check to ensure that MPI is available
    and that the number of processes is at least 2.  Assuming
    that ArrakisND is being run on a system with more than one
    CPU core.
    """
    logger = Logger("arrakis_runner", output="both")

    try:
        comm = MPI.COMM_WORLD
        size = comm.Get_size()

        # check that size >= 2, which is required
        # for Arrakis to run.  Otherwise we quit early
        if size < 2:
            logger.error(f"Number of processes must be >=2! Received {size} from MPI", "value")
    except Exception as e:
        logger.error(f"Error occurred with gathering MPI: {e}", "runtime")

    parser = argparse.ArgumentParser(
        prog="Arrakis Module Runner",
        description="This program runs the Arrakis module " + "from a config file.",
        epilog="...",
    )
    parser.add_argument(
        "config_file",
        metavar="<config_file>.yml",
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

    config_file = args.config_file
    name = args.name
    number_of_files = args.number_of_files

    try:
        config = ConfigParser(config_file).data
    except Exception as e:
        logger.error(f"failed to parse config: {e}")

    if number_of_files is not None:
        try:
            number_of_files = int(number_of_files)
        except Exception as e:
            logger.error(f"unable to convert number_of_files from {type(number_of_files)} to int: {e}")

        if isinstance(number_of_files, int):
            if number_of_files > 0:
                config["arrakis_nd"]["number_of_files"] = number_of_files
    meta = {
        "name": name,
        "config_file": config_file
    }

    try:
        arrakis = Arrakis(config, meta)
    except Exception as e:
        logger.error(f"failed to construct Arrakis object: {e}")

    arrakis.run_arrakis_nd()


if __name__ == "__main__":
    run()
