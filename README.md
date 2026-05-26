# 2D Ising Model: Phase Transitions and Universality

This repository contains the source code, data, and scripts for a computational study of the 2D Ising model. We analyze three different lattice topologies using Monte Carlo simulations (Metropolis algorithm) to study their thermodynamic properties and demonstrate the Universality hypothesis.

## 📂 Repository Structure

### Root Files
* **`plot1d.py`**: A simple script comparing the 1D vs 2D Ising model. It shows why 1D lacks a stable ferromagnetic phase at $T>0$.
* **`plot_ising.py`**: Generates an animated GIF (`ising_termometro.gif`) showing the evolution of magnetic domains as temperature changes.
* **`LICENSE`**: GNU GPLv3 License.

### `NNN/` (Square Lattice: NN & NNN)
Contains the simulation for the classic Square Lattice (Nearest Neighbors) and the extended version (Next-Nearest Neighbors).
* **`.f90` / `.f` files**: Fortran source code for the Monte Carlo simulations.

### `Triangle/` (Triangular Lattice)
Contains the adaptation of the model for a triangular topology (each spin interacts with 6 neighbors).
* **`Triangulo.f90`**: Fortran source code for the triangular lattice.
