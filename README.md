# tesi-fem

Implementation of f.e.m method in python 


## Installation

Clone the repository, create a virtual enviroment

```bash
  python3 -m venv v_env
  source v_env/bin/activate
  pip install -r requirements.txt
```


    
## Usage/Examples
```bash
   python3 main.py
```

This is a very basic version of the program, to change the problems' parameters modify the variables inside *main.py*.
When ran, the first graph shows the function **T(x)** obtained by solving the ODE with a native scipy function.
The second graph pops up once the first one is closed and it shows **T(x)** obtained by FEM method.
