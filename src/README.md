# Source code folder

All re-usable source code for the project goes here.


Different possible methods exist for checking the degree of induced irregularity - the increase in the variance of the radii of the disk.

```
number  = .no
position = .x
velocity = .v
mass = .m
forceful = .forceful
angular momentum = .AM()
kinetic energy = kin_energy()

init
X
V
AM()
kin_energy()
gravp_energy()
rk4o()
force_fun()

circular_position
circular_velocity
zmf_transform
centre_transform
tot_energies
tot_am
```

The source folder is structured as follows:


```
src
├── __init__.py    <- Makes src a Python module
│
├── constants.py   <- Includes project wide constants for easy imports
│
├── data_loading   <- Scripts to download or generate data
│
├── models         <- Scripts to train models and then use trained models to make
│                     predictions
└── tests          <- Scripts for unit tests of your functions
```

This generic folder structure is useful for most project, but feel free to adapt it to your needs.
