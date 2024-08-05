### Running the Event Display
Currently there is a script which will install ArrakisND and launch the event display from within NERSC.  After logging into NERSC, download the ArrakisND package:
```bash
git clone https://github.com/Neutron-Calibration-in-DUNE/ArrakisND 
```
Then, load python and run the script to launch the display:
```bash
module load python
cd ArrakisND
git fetch
git checkout -b develop remotes/origin/develop
source scripts/run_display.sh
```