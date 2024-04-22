### Dummy particulate simulation and sensor emulation

This is a pair of programs that generate pariculate flow data. The physical system is a source with convection-diffusion.

## Installation
`make` should work without any fuss on most linux-like environments. `getopt` is the only potentially controversial dependency.

`graph_dat.py` requires numpy and matplotlib

## dummy_sim

Simulate convection-diffusion with a source. On a user-defined interval, the solution is output to a binary file, which can be visualized with `graph_dat.py`.

Usage: `./dummy_sim -i config.toml`

## dummy_sense

Generates a time series of environmental readings of a less-deterministic flow than `dummy_sim`. The output is a series of lines, each line being a sensor reading. The sensor readings are of the format:

`{ t: <timestep>, loc: (<longitude>, <latitude>), value: <concentration>}`

Usage: `./dummy_sense -i config.toml`

## Configuration

Configuration is shared between the programs in a TOML file. Some simulation parameters can be set by option flags, but they will be overwritten by toml config, if present. Use of flags is not recommended, to the point that I am not documenting them here :-)

Example TOML file with reasonable (stable) values and annotation:

`config.toml`
```
[model]               # Relevant to both sim and sense
plume = [-39, 46.2]   # plume (Source) location [longitude, latitude]
wind = [15, 25]       # wind vector
baseline = 2          # baseline conecntration - affecs starting and boundary conditons
source = 20           # magnitude of Source term at plume coordinates
steps = 100_000       # how many steps to run
dt = 0.00001          # time interval per step
diffusivity = 1       # coefficient of diffusion

[model.grid]          # Relevant to both sim and sense
min = [-40, 45]       # lower bounds of simulation
max = [-38, 48]       # upper bounds of simulation
delta = [0.05, 0.05]  # grid deltas

[sim]                 # Relevant to sim
out_steps = 1_000     # write an output file every out_steps steps
out_dir = "sim_out"   # directory in which to place outputs

[environment]         # Relevant to sense
wind_shift = [1, 1]   # Wind will vary randomly and uniformly within (wind - wind_shift, wind + wind_shift)
variation = 0.05      # Concentration will be modified randomly within +/- variation %

[sensors]             # Relevant to sense
count = 100           # How many sensors are there
interval = 100        # How often (in steps) is sensor output generated
rate = 10             # 0-100 probability that any given sensor will generate a reading in a given interval
mobile = 10           # How many sensors are mobile
speed = 8             # Each grid square is divided into a speed x speed subgrid. A mobile sensor has a 50% chance of taking one cadrinal step on the subgrid every interal
noise = 0.05          # How noisy are the sensors readings? The reading value will vary uniformly withint +/- noise of the actual value

```


## Notes/Warnings

Absolutely no stability checking. The use of Euler's method makes it very, very easy to get instablility at otherwise reasonable model paramters. Porosity is 1.

Stability conditions can be found within this [article](https://en.wikipedia.org/wiki/Numerical_solution_of_the_convection–diffusion_equation)
