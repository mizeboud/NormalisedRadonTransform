Created by Maaike Izeboud (m.izeboud@tudelft.nl | @izeboudmaaike )

This repository provides code for the Normalised Radon Transform Damage (NeRD) detection method.
The code in this repository accompanies the paper "Damage Detection on Antarctic Ice Shelves
using the Normalised Radon Transform" by M. Izeboud and S. Lhermitte, in review at Remote Sensing of Environment (2022).

# About
NeRD is an automated method that detects damage and its orientation from a (high resolution) image.

It has so far been tested/applied on:
- Optical imagery from Sentinel-2 and LandSat 7/8 satellites
- SAR imagery from Sentinel-1 and RAMP RadarSat

NeRD follows the following steps:
- Read image (.tif) and convert to grayscale
- Cut image into small window blocks
- Apply the Normalised Radon Transform to each window
- Calculate damage signal and orientation


# Structure of the repo
This repository provides NeRD code for both Matlab and Python.

- config_files: this folder contains configuration files that are used to define parameter settings used in NeRD (python application)
- data: this folder contains an example image.tif and its detected damage
- files: this folder contains necessary files to run NeRD
- matlab: code of NeRD method, including example script
- python: code of NeRD method, including example notebook

# Data
File convention:
If an image with filename image.tif is used, NeRD output will have the following convention:
- image.nc or image.mat (python / matlab application), containing all output variables.
- image_[VARIABLE].tif (output geotiff for variables)


# How to apply NeRD
- make sure you have an image.tif available on your local machine (or use the example)
- matlab: open the script ./matlab/scripts/NeRD_example.m
    - this requires the ./matlab/functions/ directory to be part of your matlab path
    - you should update some path locations in the script
    - run the script
- python:
    - create a config.ini file in ./config_files/ with the desired settings     
    - run python from the terminal:
    `python ./scripts/main.py ./config_files/config_file.ini ./path/to/image.tif`

# Example
- ...

# Matlab code
