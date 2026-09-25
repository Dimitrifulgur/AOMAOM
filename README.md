# Angular Overlap Model (AOM) Simulator

This project is a high-performance quantum chemistry simulator written in Julia, specifically designed to calculate the multiplet energy levels of transition metal complexes. It builds the full electronic Hamiltonian in the basis of Slater Determinants and solves for the eigenstates considering interelectronic Coulomb repulsion, Spin-Orbit coupling, external magnetic fields (Zeeman effect), and Crystal Field/Angular Overlap Model (AOM) ligand perturbations.

## Features

- **Full Configuration Interaction (FCI) within the d-shell**: Generates the complete set of microstates (Slater Determinants) for a given number of $d$-electrons (e.g., $d^7$ for Co(II)).
- **Coulomb Repulsion**: Calculates precise two-electron Coulomb integrals ($1/r_{12}$) parameterized by Racah/Condon-Shortley parameters ($F_0, F_2, F_4$). Employs strict $M_L$ and $M_S$ selection rules to speed up matrix assembly.
- **Spin-Orbit Coupling**: Evaluates the $\hat{L} \cdot \hat{S}$ operator in the Configuration State Function (CSF) basis, heavily optimized via direct matrix transformations (BLAS) and $M_J$ conservation selection rules.
- **Angular Overlap Model (AOM)**: Defines arbitrary ligand geometries (e.g., octahedral, tetrahedral, trigonal prismatic) and calculates crystal field splitting based on $e_\sigma$, $e_\pi$ parameters using Wigner rotation matrices and spherical harmonics.
- **Zeeman Effect**: Includes magnetic field perturbations $\vec{B} \cdot (\hat{L} + g_e \hat{S})$ with $M_J$ selection rules for EPR/magnetism simulations.
- **Blazing Fast**: The codebase avoids standard bottleneck loops $O(N^4)$ by computing interactions in the determinant basis $O(N^2)$ and mathematically transforming them to the CSF basis via highly optimized linear algebra (`mul!`).

## Project Structure

The project has been refactored into a clean modular architecture:

- **`Main.jl`**: The entry point. It orchestrates the physics simulation. It defines the ligand geometry, builds the initial Hamiltonian, diagonalizes the Coulomb matrix, applies AOM and Spin-Orbit perturbations, and finally outputs the energy spectrum.
- **`Types.jl`**: Contains core data structures (`SpinOrbital`, `SlaterDet`, `CSF`, `Term`) mapping the mathematical representations of quantum states to Julia's type system.
- **`Integrals.jl`**: The physics engine. Computes all one- and two-electron integrals (`H_Coulomb`, `H_SO`, `H_AOM`, `H_Zeeman`) using Condon-Slater rules. Uses advanced selection rules to skip non-physical matrix elements instantly.
- **`MathUtils.jl`**: Provides algebraic utilities like `Wigner 3-j` symbol calculations, phase generation, spin matrices, and analytical $d$-orbital geometric transformations.

## Usage

1. Open `Main.jl` and configure your system parameters:
   - `nel`: Number of $d$-electrons (e.g., `nel = 7`).
   - `L`: Orbital angular momentum (e.g., `L = 2` for $d$-orbitals).
   - `F0, F2, F4`: Condon-Shortley parameters for Coulomb repulsion (in $cm^{-1}$).
   - `ligands, sigma`: Geometry definition (e.g., `get_octahedral()`).
   - Spin-Orbit coupling constant $\lambda$ (e.g., `120.0` $cm^{-1}$).

2. Run the code from the terminal:
   ```bash
   julia --project=. Main.jl
   ```

3. The script will output:
   - **Coulomb energies**: Energy levels of the free ion $LS$-terms split strictly by interelectronic repulsion.
   - **Calculated energies**: The final energy spectrum (in $cm^{-1}$) representing the fine structure of the complex after the AOM crystal field and Spin-Orbit coupling have broken the degeneracy.

## Dependencies
- `LinearAlgebra`
- `Combinatorics`
- `WignerSymbols`
- `StaticArrays`
- `SparseArrays`
- `PlotlyJS` (for optional visualization)

Ensure you have instantiated the environment using the included `Project.toml`.
