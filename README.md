# Heat and Flow Simulations

This repository contains MATLAB codes for simulating 1D heat diffusion and 2D potential flow problems, using numerical methods commonly applied in aerospace and mechanical engineering. The repository includes projects focused on heat transfer in shuttle tiles and potential flow over airfoils.

---

## Table of Contents

1. [Project Descriptions](#project-descriptions)
2. [Requirements](#requirements)
3. [Getting Started](#getting-started)
4. [Code Details](#code-details)
    - [Heat Diffusion Simulations](#heat-diffusion-simulations)
    - [Potential Flow Simulations](#potential-flow-simulations)
5. [Results and Visualizations](#results-and-visualizations)
6. [Contributing](#contributing)
7. [License](#license)

---

## Project Descriptions

### 1. Heat Diffusion Simulations
Simulates 1D heat transfer through a shuttle tile during re-entry using the heat diffusion equation. Two numerical schemes are implemented:

- **Explicit (FTCS):** Predicts bondline temperature progression and calculates time to reach 600 K.
- **Implicit:** Solves the problem with larger time steps for stability and accuracy.

Features:
- Dirichlet and Neumann boundary conditions.
- Visualizations of temperature distributions over time.
- Comparison of grid resolutions (101 and 201 points).

### 2. Potential Flow Simulations
Simulates 2D linearized potential flow over a 6% thick circular arc airfoil. Two configurations are implemented:

- **Uniform Grid (411 x 201):** Uses Gauss-Seidel relaxation to compute steady-state flow.
- **Non-Uniform Grid (141 x 51):** Refines the grid spacing near the airfoil for increased accuracy.

Features:
- Steady-state solutions for elliptic PDEs.
- Accurate boundary conditions for symmetry planes and freestream flow.
- Contour plots of Mach number and pressure coefficient distributions.

---

## Requirements

- **MATLAB** (R2018b or later recommended)
- Basic understanding of numerical methods and fluid/thermal dynamics.

---

## Getting Started

1. Clone the repository:
   ```bash
   git clone https://github.com/username/HeatAndFlowSimulations.git
   cd HeatAndFlowSimulations
2. Open MATLAB and navigate to the repo folders.
3. Run the desired simulation:
   - For heat diffusion, open and execute `heat_diffusion_explicit.m` or `heat_diffusion_implicit.m`.
   - For potential flow, open and execute `potential_flow_uniform.m` or `potential_flow_nonuniform.m`.

4. Adjust parameters as needed:
   - Grid sizes, time steps, or material properties can be modified within the scripts to explore different scenarios.

---

## Code Details

### Heat Diffusion Simulations

#### Inputs:
- Tile thickness: 1 inch.
- Material properties: Density (140 kg/m³), Specific Heat (628 J/kg·K), Thermal Conductivity (0.048 W/m·K).
- Boundary conditions: 2000 K on the re-entry side (Dirichlet) and adiabatic (Neumann) on the bondline.

#### Outputs:
- Time to reach 600 K at the bondline (with linear interpolation for precise results).
- Temperature distributions at multiple time steps.

#### Numerical Schemes:
- **Explicit Scheme (FTCS):** A first-order temporal, second-order spatial method. Requires small time steps for stability.
- **Implicit Scheme:** A more stable method allowing larger time steps while maintaining accuracy.

#### Key Functions:
- `plotTemperatureProfiles`: Visualizes the temperature distribution at specific time steps.

---

### Potential Flow Simulations

#### Inputs:
- Airfoil: 1-meter chord length, 6% thickness, at 0° angle of attack.
- Freestream Conditions:
  - Mach Number: 0.5
  - Temperature: 300 K
  - Pressure: 1 bar (100 kPa)

#### Outputs:
- Normalized residual history for convergence.
- Contour plots of Mach number distributions.
- Pressure coefficient distributions along the airfoil.

#### Numerical Scheme:
- **Gauss-Seidel Relaxation with Acceleration:** Solves elliptic PDEs for steady-state solutions.

#### Grid Configurations:
- **Uniform Grid:** 411 x 201 grid points with fixed spacing for simplicity.
- **Non-Uniform Grid:** 141 x 51 grid points with refined spacing near the airfoil for improved accuracy.

---

## Results and Visualizations

### Heat Diffusion
- **Temperature vs. Time:** Shows the progression of the bondline temperature reaching 600 K.
- **Temperature Profiles:** Multiple plots showing the distribution of temperature across the tile at various time steps.

### Potential Flow
- **Convergence Plot:** Residual history normalized to the initial iteration to demonstrate solution stability.
- **Mach Number Contour:** Visualizes flow patterns around the airfoil.
- **Pressure Coefficient Plot:** Displays the pressure distribution along the airfoil surface.

---

## Contributing
This was for Dr. Bill Engblom AE 455 class, as such the codes are based on the codes he provided in class.

---


