from setuptools import find_packages
from setuptools import setup

with open("README.md", "r") as file:
    long_description = file.read()

setup(
    # name
    name="arrakis_nd",
    # current version
    #   MAJOR VERSION:  02
    #   MINOR VERSION:  00
    #   Maintenance:    00
    version="2.00.00",
    # descriptions
    description="Arrakis module for near detector data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords="",
    # my info
    author="Nicholas Carrara, Marjolein van Nuland, Luis Lepin",
    author_email="ncarrara.physics@gmail.com",
    # where to find the source
    url="https://github.com/Neutron-Calibration-in-DUNE/ArrakisND",
    # requirements
    install_reqs=[],
    # packages
    packages=find_packages(
        # where='pointnet',
        exclude=["test"],
    ),
    include_package_data=True,
    # classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Experimental Physics",
        "License :: GNU",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">3.7",
    # possible entry point
    entry_points={
        "console_scripts": [
            "arrakis_nd = arrakis_nd.programs.run_arrakis:run",
            "create_arrakis_runs = arrakis_nd.programs.create_arrakis_runs:run",
            'create_plugin = arrakis_nd.programs.plugin_creator:run',
        ],
    },
)
