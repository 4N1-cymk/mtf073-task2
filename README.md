# CFD Task 2 – 2D Convection–Diffusion Solver (Python)

This repository contains my solution for **Task 2** of the Computational Fluid Dynamics course.  
The project implements a **2D convection–diffusion temperature solver** on a **non-uniform structured mesh** using the **finite volume method (FVM)**.

---

## 📌 Features

- 2D steady & unsteady temperature solver  
- Non-equidistant structured mesh  
- Hybrid convection scheme  
- Implicit time discretization  
- TDMA (line-by-line) & Gauss–Seidel solvers  
- Automatic plotting:
  - Temperature contours  
  - Heat flux vectors  
  - Residual history  
  - Mesh and velocity field  
  - Probe-point time evolution (unsteady)  
  - Optional animated GIF of the temperature field  

---

## 📂 Repository Structure

```
.
├── T2_template.py                     # Main driver script
├── T2_codeFunctions_template.py       # Numerical methods and FVM functions
├── T2_test_codeFunctions_template.py  # Unit tests for each function
├── meshAndVelocityData/               # Mesh + velocity .npz files
├── refData/                           # Reference data for testing
├── Figures/                           # Auto-generated figures
└── README.md
```

---

## ▶️ How to Run

## ⚙️ Main Settings (inside `T2_template.py`)

```python
caseID = 24            # Select case
grid_type = 'coarse'   # or 'fine'
unsteady = True        # True = transient, False = steady

deltaT = 1.0
endTime = 400.0

solver = 'TDMA'        # 'GS' or 'TDMA'
createAnimatedPlots = True
```

You can modify:
- Boundary conditions  
- Material properties  
- Probe points  
- Grid type (coarse/fine)  

---

## 📊 Output

The script automatically generates:

- Temperature contour plot  
- Mesh plot  
- Velocity vector field  
- Heat flux vectors  
- Residual convergence plot  
- Time evolution at probe points (unsteady)  
- Optional `animated_contour.gif` showing temperature evolution  

All outputs are saved inside the `Figures/` folder.

---

## 📝 Notes

- All FVM and solver logic is implemented inside `T2_codeFunctions_template.py`.
- Mesh and velocity fields are loaded from `.npz` files in `meshAndVelocityData/`.
- The test script compares your computed arrays to reference arrays.
