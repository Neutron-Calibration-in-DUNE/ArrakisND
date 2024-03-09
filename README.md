# ArrakisND

[![ArrakisCI](https://github.com/BlipDev/ArrakisCI/actions/workflows/arrakis-ci.yml/badge.svg?branch=master)](https://github.com/BlipDev/ArrakisCI/actions/workflows/arrakis-ci.yml)
[![CodeFactor](https://www.codefactor.io/repository/github/neutron-calibration-in-dune/arrakisnd/badge)](https://www.codefactor.io/repository/github/neutron-calibration-in-dune/arrakisnd)
[![Flake8 Status](./.reports/flake8/flake8-badge.svg?dummy=8484744)](./.reports/flake8/index.html)
[![Coverage Status](./.reports/coverage/coverage-badge.svg?dummy=8484744)](./.reports/coverage/index.html)

Arrakis module for the near detector data. This module generates ARRAKIS.h5 files which correspond to FLOW.h5 files whose information can be found on the [2x2 Wiki](https://github.com/DUNE/2x2_sim/wiki).  Arrakis is built to run with h5 + MPI through the use of a [dedicated container]() which is available on NERSC.  Arrakis operates as a post FLOW process which generates (1) arrays that compliment the reconstructed charge and light datasets in the FLOW, as well as (2) a set of StandardRecord objects that can be used as a benchmark for any downstream reconstruction task.

#### Developers
Nicholas Carrara, UC Davis, Dept. of Physics [nmcarrara@ucdavis.edu],

Marjolein van Nuland, National Institute for Subatomic Physics (NIKHEF)[mnuland@nikhef.nl]

### Table of Contents

1. [ Getting the Repository ](#get)
2. [ Getting the Container ](#container)
3. [ Running on Perlmutter ](#perlmutter)
4. [ Quick Start ](#quickstart)
5. [ Labeling Logic ](#labelinglogic)
6. [ Usage ](#usage)
7. [ Versioning ](#versions)
8. [ Contact (Authors) ](#contact)
9. [ Citation ](#citation)
10. [ License ](#license)
11. [ Support ](#support)

<a name="get"></a>
## Getting the Repository

In the terminal, one can clone this repository by typing the command:

`git clone https://personal_username@github.com/Neutron-Calibration-in-DUNE/ArrakisND.git`

This uses the HTTPS protocol. For environments (e.g. computing clusters) where one has to use the SSH protocol:

`git clone git@github.com:Neutron-Calibration-in-DUNE/ArrakisND.git`

Anyone in the "Neutron-Calibration-in-DUNE" organization should be able to develop (push changes to the remote repository).

Please contact Nicholas Carrara, Marjolein van Nuland or Georgette Kufatty about becoming involved in development before merging with the master branch. 

<a name="container"></a>
## Getting the Container
Arrakis is developed to work within a container environment.  The container is hosted on [dockerhub](https://hub.docker.com/r/infophysics/nersc), which can be downloaded through the command:

```bash
docker pull infophysics/nersc
```

<a name="perlmutter"></a>
## Running on Perlmutter
Running on Perlmutter requires that you have access to NERSC and to the 'dune' group.  

### Setting up a daily SSH
Below is a snippet taken from the [NERSC site](https://docs.nersc.gov/connect/mfa/#sshproxy) which details how to set up a daily ssh proxy, so that you only have to enter your password + OTP once per day to access perlmutter:

NERSC has developed a service, called sshproxy, that allows you to use MFA to get an ssh key that is valid for a limited time (24 hours by default). sshproxy provides a type of single-sign-on capability for ssh to NERSC systems. Once you have obtained a key, you can use it to ssh to NERSC systems without further authentication until the key expires.

The sshproxy service uses a RESTful API for requesting keys. NERSC provides a bash client script that you can use from the command line on a Unix-like computer.

sshproxy on Unix-like Systems (macOS, Cygwin and Windows Subsystem for Linux included)¶
INSTALLING THE CLIENT¶
You can download the bash client sshproxy.sh via scp:

```bash
scp myusername@dtn01.nersc.gov:/global/cfs/cdirs/mfa/NERSC-MFA/sshproxy.sh .
```

where myusername is your NERSC login ID. The above command uses a data transfer node (dtn01), but you can use any machine which you can access that can access the Community file system.

USING SSHPROXY¶
The sshproxy client, without any arguments, will use your local username, and obtain an ssh key with the default lifetime (24 hours). The private and public key will have the names nersc and nersc-cert.pub, and will be stored in your ~/.ssh directory.

Run the sshproxy.sh script from where you installed it. The script will prompt you to enter your password and OTP, in the same manner as you would do to ssh to a NERSC system with MFA:
```bash
$ ./sshproxy.sh -u <nersc_username>
```
Enter your password+OTP:
Enter your NERSC password immediately followed by OTP as a single string, as before. Upon successfully authenticating, the client will install an ssh key and display a message showing the path to the key pair installed on your local computer and the expiration date and time for the keys. By default, the name of the files will be ~/.ssh/nersc and ~/.ssh/nersc-cert.pub (you can change the name with a command-line argument).

### Loading the Container



<a name="quickstart"></a>
## Quick Start

<a name="labelinglogic"></a>
## Labeling Logic
ArrakisND works in the same spirit as Arrakis, which is the LArSoft version of this software.    

| Feature | Type | Description |
| ------- | ---- | ----------- |
| x | float | position of a charge signal along the x direction | 
| y | float | position of a charge signal along the y direction |
| z | float | position of a charge signal along the z direction |
| Q | float | value of the reconstructed charge signal |

The second stage consists of a set of functions, each of which processes different particles and their detector output.  The labeling scheme currently consists of assigning seven different labels to each reconstructed charge/light point in the TPC.  The seven labels are shown in the following tables:

| High-level features | Description |
|---------------------|-------------|
| topology            | topological descriptor of physics (e.g. blip, track, shower) |
| particle            | the pdg code of the particle which caused the energy deposition |
| physics_micro       | low-level descriptor of physics processes (e.g. mip_ionization, gamma_conversion, etc.) |
| physics_meso        | mid-level descriptor of physics processes (e.g. hip, delta_electron, capture_gamma, etc.) |
| physics_macro       | high-level descriptor of physics processes (e.g. cc_qe, neutron_capture, radiological, etc.) |
| unique_topology     | unique identifier of individual topology instances |
| unique_particle     | unique identifier of individual particle instances (i.e. track id) |
| unique_physics_*    | unique identifier of individual physics instances |

For information on how each of these labels is assigned to each reconstructed point, see the documentation for [BLIP](https://github.com/Neutron-Calibration-in-DUNE/Blip) which can be found at [![Documentation Status](https://readthedocs.org/projects/blip-dune/badge/?version=latest)](https://blip-dune.readthedocs.io/en/latest/?badge=latest)

<a name="usage"></a>
## Usage

<a name="versions"></a>
## Versioning
For the versions available, see the [tags on this repository](https://github.com/Neutron-Calibration-in-DUNE/ArrakisND/tags). 
   
<a name="contact"></a>
## Contact (Authors)
If you have questions, please contact 

Nicholas Carrara, nmcarrara@ucdavis.edu,

Marjolein van Nuland, mnuland@nikhef.nl

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
