"""
ArrakisND main program
"""
import argparse

from arrakis_nd.utils.logger import Logger
from arrakis_nd.utils.display.arrakis_display import ArrakisDisplay


def run():
    """
    This program runs the Arrakis display and generates a 
    dash server instance that can be run on NERSC ...
    """

    """Set up command line arguments"""
    parser = argparse.ArgumentParser(
        prog="Arrakis Display Runner",
        description="This program runs the Arrakis display and generates a \
            dash server instance that can be run on NERSC ... ",
        epilog="...",
    )

    """Parse command line arguments"""
    args = parser.parse_args()
    
    """Create the Arrakis instance"""
    try:
        arrakis = ArrakisDisplay()
    except Exception as e:
        raise RuntimeError(f"failed to construct Arrakis object: {e}")

    """Run Arrakis"""
    arrakis.run_app()


if __name__ == "__main__":
    run()
