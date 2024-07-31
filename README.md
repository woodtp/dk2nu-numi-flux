# Simple Dk2Nu Reader
The purpose of this project is to read [G4NuMI](https://cdcvs.fnal.gov/redmine/projects/numi-beam-sim/wiki) flux files in the Dk2Nu format, calculate the weights for the input position, and write a tree to a new ROOT file.
Output branches are defined in `spectra_definitions.py`, and the following are the default branches in the output file:
```
nu_pdg                     -- Neutrino PDG code
parent_pdg                 -- Neutrino's Parent PDG code
weight                     -- Flux weight for the input position
nu_energy                  -- Neutrino energy [GeV]
vx                         -- Neutrino production vertex x-position in NuMI coordinates [cm]
vy                         -- Neutrino production vertex y-position in NuMI coordinates [cm]
vz                         -- Neutrino production vertex z-position in NuMI coordinates [cm]
ppvx                       -- Neutrino's parent production vertex x-position in NuMI coordinates [cm]
ppvy                       -- Neutrino's parent production vertex y-position in NuMI coordinates [cm]
ppvz                       -- Neutrino's parent production vertex z-position in NuMI coordinates [cm]
pdpx                       -- Neutrino's parent momentum at the decay point; x-component [GeV/c]
pdpy                       -- Neutrino's parent momentum at the decay point; y-component [GeV/c]
pdpz                       -- Neutrino's parent momentum at the decay point; z-component [GeV/c]
parent_momentum            -- Neutrino's parent momentum at the decay point [GeV/c]
theta_p                    -- Neutrino's parent momentum angle with respect to the beam direction, in the plane formed by the beam and the input position [rad]
ancestor_parent_pdg        -- Vector of the PDG codes for each of the neutrino's ancestors. The first element is the beam proton.
ancestor_parent_mom        -- Vector of momenta for each of the neutrino's ancestors. The first element is the beam proton; excludes the neutrino. [GeV/c]
ancestor_produced_mom      -- Vector of momenta for each of the neutrino's ancestors, shifted by 1 position. The first element is the particle produced by the beam proton's initial interaction; excludes the neutrino. [GeV/c]
ancestor_produced_theta    -- Vector of angles between each of the neutrino's ancestors momenta and the corresponding incident particle. The first element is the particle produced by the beam proton's initial interaction; excludes the neutrino. [rad]
ancestor_pT                -- Vector of the transverse component of momenta for each ancestor in frame of the interaction which produced it. [GeV/c]
ancestor_xF                -- Vector of the Feynman-x for each ancestor in frame of the interaction which produced it.
ancestor_vol               -- Vector of the interaction volumes for each ancestor. The first element is the volume of the beam proton's initial interaction; usually `TGT1` corresponding to the carbon target.
```

## Instructions to build

### OS-Independent Installation
If doing an OS-Independent installation, the following dependencies are required:

Prerequisites:
- [ROOT](https://root.cern/install/) >=v6
- [VDT](https://github.com/dpiparo/vdt) >=v04.4

Both are prerequisites for [Dk2Nu](https://github.com/NuSoftHEP/dk2nu), which will built as part of the project.

```shell
$ git clone https://github.com/woodtp/dk2nu-numi-flux.git  # Clone the repository
$ cd dk2nu-numi-flux  # Move into the project directory
$ ./build.sh  # Build the project
```

### Building in an ICARUS GPVM
The `build_sl7.sh` is intended to be *sourced* Scientific Linux Environment in an ICARUS GPVM.
This script will set up the following dependencies:
- CMake v3_27_4
- Dk2Nu v01_10_01
- ROOT v6_28_10a

```shell
$ source ./build_sl7_sourceme.sh  # Build the project, have to source instead of execute if on a FNAL GPVM
```

## Running the code
```shell
# If on ICARUS GPVM, source the sl7 setup script
$ source setup_sl7_sourceme.sh

# Otherwise, if using the OS-independent installation, then just activate the python environment
$ source venv/bin/activate

$ ./flux_ana.py  # Run the code
```

### Usage
```
usage: flux_ana.py [-h] [-c CONFIG] [-f] [-d]

Reads Dk2Nu format, calculates weights for the input position, and writes a
tree to a new ROOT file.

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
                        Path to configuration file. If unspecified, defaults to the `config.toml` in the current directory.
  -f, --overwrite       Overwrite output file if it exists
  --mt                  Use multithreading
  -d, --debug           run in debug mode
```

## Configuration File

```toml
# Sample config.toml

location = [450.37, 7991.98, 79512.66]  # ICARUS TPC center in NuMI coords (x, y, z) [cm]

output_file = "output.root"             # Name of the output file

pot_per_file = 500_000

[file_sets]
nominal = "/path/to/files/*.root"   # Path to flux files in Dk2Nu format
# NOTE: The name of the file set will be used as the key for creating TTrees in the output file
```

# References
- [G4NuMI](https://cdcvs.fnal.gov/redmine/projects/numi-beam-sim/wiki)
- [Dk2Nu](https://github.com/NuSoftHEP/dk2nu)
- [ROOT](https://root.cern/install/) >=v6
- [VDT](https://github.com/dpiparo/vdt) >=v04.4
