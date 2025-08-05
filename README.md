# A Python code for solving the Helmholtz equation in a two-dimensional waveguide.

We use NGSolve (https://ngsolve.org/) to solve the Helmholtz equation in a two-dimensional waveguide.  In the provided example, we solve the Helmholtz equation in a waveguide with a variable depth, and specifically a waveguide of increasing depth. A single circular scatterer is placed in the waveguide.

## Installation

### Required Packages

The following packages are required to run this code:
- `numpy` - Numerical computing library
- `scipy` - Scientific computing library  
- `meshio` - I/O for various mesh formats

Install the required packages:

```bash
pip install numpy scipy meshio
```

### NGSolve Installation

To install NGSolve, run the following command:
```bash
pip install --upgrade ngsolve
```
For alternative installation methods, see https://ngsolve.org/downloads.

## Usage

```bash
netgen HE_in_WG_var_depth.py
```
or

```bash
python3 HE_in_WG_var_depth.py
```
Running the code will generate the solution and plot it. It also saves two files inside the `data` folder:
- `mesh_wg_var_depth.msh`: the mesh file.
- `sol_wg_var_depth_f73.0.mat`: the solution file.

Both files can be used to visualize the solution outside of NGSolve.

