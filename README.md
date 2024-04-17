# NuMI Bug Studies

## Instructions to build
```sh
$ . ./build.sh  # Build the project, have to source instead of execute if on a FNAL GPVM
```

## Running the code
```sh
$ . setup.sh  # Setup the environment, have to source instead of execute if on a FNAL GPVM
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
                        path to configuration file
  -f, --overwrite       overwrite output file if it exists
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
