from setuptools import find_packages
from setuptools import setup

with open("README.md", "r") as file:
    long_description = file.read()

# Function to read the list of requirements from requirements.txt
def read_requirements():
    with open('requirements.txt') as req:
        return req.read().splitlines()

major_version = 2
minor_version = 0
maintenance = 2

setup(
    # name
    name="arrakis_nd",
    # current version
    version=f"{major_version}.{minor_version}.{maintenance}",
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
    install_requires=read_requirements(),
    # packages
    packages=find_packages(
        exclude=["test"],
    ),
    include_package_data=True,
    # classifiers
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">3.7",
    # possible entry point
    entry_points={
        "console_scripts": [
            "arrakis_nd = arrakis_nd.programs.run_arrakis:run",
            "arrakis_display = arrakis_nd.programs.run_arrakis_display:run",
            "create_arrakis_runs = arrakis_nd.programs.create_arrakis_runs:run",
            'create_plugin = arrakis_nd.programs.plugin_creator:run',
        ],
    },
)
