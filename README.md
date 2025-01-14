# Atomic Lattice Simulation for Sound Speed in Copper

## Overview
This project simulates the motion of atoms in a one-dimensional copper lattice to calculate the speed of sound in copper. The simulation leverages physical constants and parameters such as atomic mass, density, and Young's modulus to accurately model the system.

## Features
- **Atomic Mass Calculation**: Utilizes the atomic mass unit to determine the mass of a copper atom.
- **Lattice Initialization**: Configures a lattice of atoms with initial positions, velocities, and accelerations.
- **Harmonic Oscillator Dynamics**: Implements Euler's method to update atomic positions over time.
- **Damping Coefficient**: Incorporates a damping factor to simulate realistic atomic interactions.
- **Visualization**: Generates and saves images of the atomic lattice at various time steps.
- **Sound Speed Calculation**: Records and plots the calculated speed of sound over time.

## Simulation Parameters
- **Number of Atoms (N)**: 40
- **Lattice Constant (L0)**: Derived from the atomic mass and density.
- **Spring Constant (k)**: Calculated using Young's modulus and the lattice constant.
- **Damping Coefficient (c)**: Set to ensure realistic damping of atomic motion.
- **Time Step (dt)**: Chosen to be small relative to the harmonic oscillator period.
- **Initial Displacement (dx)**: Increased to ensure noticeable atomic movement.

## Usage
The simulation runs until the last atom in the lattice has moved 1% from its equilibrium position, indicating that the wave has propagated through the entire lattice.

## Output
- **Sound Speed Data**: Saved to a file for further analysis.
- **Visualization Images**: Generated at specified intervals to visualize atomic motion.

## Dependencies
- **NumPy**: For numerical computations.
- **Matplotlib**: For plotting and visualization.

## Installation
To install the required dependencies, run:
```bash
pip install numpy matplotlib
