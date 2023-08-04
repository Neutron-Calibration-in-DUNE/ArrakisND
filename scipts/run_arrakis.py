"""
"""
import argparse

from arrakis_nd.arrakis import Arrakis

parser = argparse.ArgumentParser(
    prog='ArrakisND',
    description='This program constructs a ArrakisND module ' +
        'from a config file, and then runs training/inference.',
    epilog='...'
)
parser.add_argument(
    'config_file', metavar='<str>.yml', type=str,
    help='config file specification for a ArrakisND module.'
)
parser.add_argument(
    '-n', dest='name', default='arrakis',
    help='name for this run (default "arrakis").'
)

if __name__ == "__main__":
    
    args = parser.parse_args()
    arrakis_module = Arrakis(
        args.config_file
    )