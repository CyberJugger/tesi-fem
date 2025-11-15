# tesi-fem

Implementation of F.E.M method in Python.

## Installation

Clone the repository and create a virtual environment:

```bash
python3 -m venv v_env
source v_env/bin/activate
pip install -r requirements.txt
# Linux: sudo apt-get install python3-tk
# macOS: brew install python-tk
```

## Usage / Examples

The repository contains two main folders:

```bash
calore\
timoshenko\
```

### Heat Problem

```bash
python3 calore/gui.py
```

The interface shows several fields, organized into three sections:

**Domain**

* **x_start** → starting coordinate of the domain 
* **x_end** → ending coordinate of the domain
* **n elements** → number of finite elements 

**Boundary Conditions**

* **left / right** → choose type (Dirichlet or Neumann) and set the corresponding value

  * **Dirichlet** → specify temperature at the boundary
  * **Neumann** → specify heat flux at the boundary

**Functions**

* **k(x)** → thermal conductivity (W/m·K), can be constant or function of `x`
* **a(x)** → cross-sectional area (m²), can be constant or function of `x`
* **f(x)** → heat source (W/m³), can be constant or function of `x`

Pressing **Run** computes the temperature distribution along the domain and visualizes it using `matplotlib`.

### Timoshenko FEM GUI

```bash
python3 timoshenko/gui_timoshenko.py
```

**Domain**

This section allows the definition of the **beam domain and mesh**:

* **x start** → starting coordinate of the beam 
* **x end** → ending coordinate of the beam 
* **n elements** → number of finite elements 

**Material / Section**

Define the beam material and cross-sectional properties:

* **E(x)** → Young's modulus (Pa), constant or function of `x`
* **I(x)** → section moment of inertia (m⁴), constant or function of `x`
* **G(x)** → shear modulus (Pa), constant or function of `x`
* **A(x)** → cross-sectional area (m²), constant or function of `x`
* **kappa** → Timoshenko shear coefficient (default `5/6`)

**Load**

* **q(x)** → distributed load along the beam (N/m), constant or function of `x`

**Boundary Conditions**

* **left w / phi** → transverse displacement and rotation at the left end; leave empty for unconstrained (None)
* **right w / phi** → transverse displacement and rotation at the right end; leave empty for unconstrained (None)

⚠️ All boundary constraints are currently **Dirichlet**; Neumann conditions (forces/moments) are not included to ensure system stability.

**Execution and Results Visualization**

Pressing **Run** calculates the beam transverse displacement `w(x)` and rotation `φ(x)`, displaying two separate plots using `matplotlib`.

**Additional Notes**

* All parameters (E, I, G, A, q) can be constants or functions of `x`, allowing modeling of beams with variable material or cross-section.
* The solver uses the 1D Timoshenko beam formulation, supporting both uniform and non-uniform beams.

