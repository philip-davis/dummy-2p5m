### Dummy particulate simulation and sensor emulation

This is a pair of programs that generate pariculate flow data. The physical system is a source with convection-diffusion.

## Installation
`make` should work without any fuss on most linux-like environments. `getopt` is the only potentially controversial dependency.

`graph_dat.py` requires numpy and matplotlib

## dummy_sim

Simulate convection-diffusion with a source. On a user-defined interval, the solution is output to a binary file, which can be visualized with `graph_dat.py`.

## dummy_sense

Generates a time series of environmental readings of a less-deterministic flow than `dummy_sim`.

## Configuration

## Notes/Warnings

Absolutely no stability checking. The use of Euler's method makes it very, very easy to get instablility at otherwise reasonable model paramters.
