# tesi-fem

Implementation of f.e.m method in python 


## Installation

Clone the repository and create a virtual environment

```bash
  python3 -m venv v_env
  source v_env/bin/activate
  pip install -r requirements.txt
  linux: sudo apt-get install python3-tk 
  macos: brew install python-tk
```


    
## Usage/Examples
```bash
   python3 gui.py
```

The interface shows various fields, including three sections:

**Dominio**
- x_start: start of the domain
- x_end: end of the domain
- n: number of elements in which the domain is devided
**Condizioni al contorno**

it has two fields allowing the user to choose between dirichlet or neumann boundary conditions and set the values

**Funzioni**
- k(x): conductivity
- a(x): for area
- f(x): function for heat source

