# Steady-State 2D Velocity Solver (Meshfree FPM) - C++11

This project implements a 2D **steady-state velocity field solver** using the **Finite Point Method (FPM)** with a Gaussian weight kernel and a second-order polynomial basis. The code is written in **C++11** and uses **Eigen** for sparse matrix operations.

---

## 🧠 Method Overview

- **Domain**: Rectangular (2.0 × 1.0)
- **Discretization**: Uniform point cloud (meshfree)
- **Neighbor search**: Voxel-based
- **Differential operator**: Weighted Least Squares (second-order polynomial)
- **Solver**: Eigen's SparseLU for solving the linear system
- **Boundary Conditions**: 
  - Dirichlet on top and bottom boundaries
  - Interior treated via Laplacian constraints
- **Output**: CSV file containing `(x, y, velocity)`

---

## 📁 Files

| File            | Description                                      |
|-----------------|--------------------------------------------------|
| `steadystate.cpp` | Main source code for the solver                |
| `task.json`       | VS Code build configuration for `Ctrl+Shift+B` |
| `Velocity.csv`    | Output file containing computed velocity field |

---

## 🧰 Requirements

- C++11 or above
- [Eigen library](https://eigen.tuxfamily.org/) (Header-only, no installation needed)
- A C++ compiler (e.g., g++, clang++)
- Visual Studio Code (recommended for easy build via `task.json`)

---

## ⚙️ How to Compile & Run

### Using VS Code:
1. Open this folder in VS Code.
2. Ensure `Eigen` is included properly (it's assumed to be in your include path).
3. Press **`Ctrl + Shift + B`** to build and run (uses `task.json`).

### Or, Compile Manually:
```bash
g++ -std=c++11 steadystate.cpp -I /path/to/eigen -o steadystate
./steadystate
