Created by Maaike Izeboud (m.izeboud@tudelft.nl | @izeboudmaaike )

This repository provides code for the Normalised Radon Transform Damage (NeRD) detection method.
The code in this repository accompanies the paper "Damage Detection on Antarctic Ice Shelves
using the Normalised Radon Transform" by M. Izeboud and S. Lhermitte, in review at Remote Sensing of Environment (2022).

Provided are functions in both Matlab and Python.

# About
NeRD is an automated method that detects damage and its orientation from a (high resolution) image.

NeRD follows the following steps:
- Read image and convert to grayscale
- Cut image into small window blocks
- Apply the Normalised Radon Transform to each window
- Calculat crevasse signal and orientation

It has so far been applied on:
- Optical imagery from Sentinel-2 and LandSat 7/8
- Radar imagery from Sentinel-1 and RAMP RadarSat
