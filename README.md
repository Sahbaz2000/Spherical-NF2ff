# Spherical-NF2FF

A complete and modular MATLAB implementation for **spherical near-field to far-field transformation** (NF2FF) using simulated near-field data.

## 📌 Overview

This project provides a functional workflow to simulate, transform, and visualize the electromagnetic radiation characteristics of antennas using spherical scanning methods. Based on the theoretical foundation from J.E. Hansen’s *Spherical Near-Field Antenna Measurements*, the code includes all steps from near-field data generation to far-field validation.

## 🛠️ Features

- ✅ Near-field simulation using MATLAB Antenna Toolbox
- ✅ Spherical NF2FF transformation (spectral method)
- ✅ Visual output: 3D diagrams, polar plots, theta/phi cuts
- ✅ Automated validation against analytically simulated far-field data
- ✅ Fully modular MATLAB functions with documentation

## 📂 Structure

📁 misc_functions/              # Core mathematical and helper functions
📁 transformation_functions/    # Main NF2FF transformation
📁 plot_functions/              # Plotting and visualization
📁 generator_functions/         # Near-field and far-field simulation
📄 nf2ff_spherical.m            # Main script (start here)
