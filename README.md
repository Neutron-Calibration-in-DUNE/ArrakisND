# ArrakisND

[![Pytest](https://github.com/Neutron-Calibration-in-DUNE/ArrakisND/actions/workflows/test_package.yml/badge.svg?branch=master)](https://github.com/Neutron-Calibration-in-DUNE/ArrakisND/actions/workflows/test.yml)

Arrakis module for the near detector data. This module can run in two modes, stand-alone or as part of an H5Flow. Input files should be flow files from nd-larsim or data runs, which are then converted into training data sets for Blip. 

#### Developers
Nicholas Carrara, UC Davis, Dept. of Physics [nmcarrara@ucdavis.edu],

Marjolein van Nuland, National Institute for Subatomic Physics (NIKHEF)[mnuland@nikhef.nl]

### Table of Contents

1. [ Getting the Repository ](#get)
2. [ Quick Start ](#quickstart)
3. [ Labeling Logic ](#labelinglogic)
4. [ Usage ](#usage)
5. [ Versioning ](#versions)
6. [ Contact (Authors) ](#contact)
7. [ Citation ](#citation)
8. [ License ](#license)
9. [ Support ](#support)

<a name="get"></a>
## Getting the Repository

In the terminal, one can clone this repository by typing the command:

`git clone https://personal_username@github.com/Neutron-Calibration-in-DUNE/ArrakisND.git`

This uses the HTTPS protocol. For environments (e.g. computing clusters) where one has to use the SSH protocol:

`git clone git@github.com:Neutron-Calibration-in-DUNE/ArrakisND.git`

Anyone in the "Neutron-Calibration-in-DUNE" organization should be able to develop (push changes to the remote repository).

Please contact Nicholas Carrara, Marjolein van Nuland or Georgette Kufatty about becoming involved in development before merging with the master branch. 

<a name="quickstart"></a>
## Quick Start

<a name="labelinglogic"></a>
## Labeling Logic
ArrakisND works in the same spirit as Arrakis, which is the LArSoft version of this software.  The basic idea is to split up the labeling into two stages, the first of which gathers all relevant event information into various maps and arrays that can be queried easily.  The TPC data is organized in the following arrays:

| Feature | Type | Description |
| ------- | ---- | ----------- |
| x | float | position of a charge signal along the x direction | 
| y | float | position of a charge signal along the y direction |
| z | float | position of a charge signal along the z direction |
| Q | float | value of the reconstructed charge signal |

The second stage consists of a set of functions, each of which processes different particles and their detector output.  The labeling scheme currently consists of assigning seven different labels to each reconstructed charge/light point in the TPC.  The seven labels are shown in the following tables:

| High-level features | Type | Description|
| ------- | ---- | ----------- |
| source | class | denotes the source of the primary interaction (e.g. beam, radiological, pns, etc.) |
| topology | class | topological descriptor of physics (e.g. blip, track, shower) |
| particle | class | the pdg code of the particle which caused the energy deposition |
| physics | class | high-level descriptor of physics processes (e.g. hip, mip, capture gamma, etc.) |
| unique_topology | cluster | unique identifier of individual topology instances |
| unique_particle | cluster | unique identifier of individual particle instances (i.e. track id) |
| unique_physics | cluster | unique identifier of individual physics instances |

For information on how each of these labels is assigned to each reconstructed point, see the documentation for [BLIP](https://github.com/Neutron-Calibration-in-DUNE/Blip) which can be found at [![Documentation Status](https://readthedocs.org/projects/blip-dune/badge/?version=latest)](https://blip-dune.readthedocs.io/en/latest/?badge=latest)

<a name="usage"></a>
## Usage

<a name="versions"></a>
## Versioning
For the versions available, see the [tags on this repository](https://github.com/Neutron-Calibration-in-DUNE/ArrakisND/tags). 
   
<a name="contact"></a>
## Contact (Authors)
If you have questions, please contact Nicholas Carrara, nmcarrara@ucdavis.edu.

See also the list of [contributors](https://github.com/orgs/Neutron-Calibration-in-DUNE/people) who participate in this project.

See AUTHORS.md for information on the developers.

<a name="support"></a>
## Support

* Bugs: Please report bugs to the [issue tracker on Github](https://github.com/Neutron-Calibration-in-DUNE/ArrakisND/issues) such that we can keep track of them and eventually fix them.  Please explain how to reproduce the issue (including code) and which system you are running on.
* Help: Help can be provided also via the issue tracker by tagging your issue with 'question'
* Contributing:  Please fork this repository then make a pull request.  In this pull request, explain the details of your change and include tests.
   
<a name="citation"></a>
## Citation

When you use `arrakis_nd`, please say so in your slides or publications (for publications, see Zenodo link above).  This is important for us being able to get funding to support this project.
