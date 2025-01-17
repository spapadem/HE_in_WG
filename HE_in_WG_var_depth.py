# Preamble, loading necessary packages
from ngsolve import *
from netgen.occ import *
from netgen.geom2d import SplineGeometry
from scipy.io import savemat
from scipy.interpolate import griddata
import numpy as np
import meshio
import ngsolve.internal as ngsint
ngsint.viewoptions.drawoutline=0 # disable triangle outline when plotting.
import netgen.gui
import numpy as np
from Gfu_to_mesh import ConvertSolutiononMesh



## Problem setup.
# Wavenumber & source.
c0 = 1500. # Constant wave speed.
f0 = 75.   # Reference frequency.
lambda_0 = c0/f0; # Reference wavelength. We use this to define all sizes with respect to this scale.


# Waveguide characteristics
Dc = 7*lambda_0 # Waveguide constant depth
Dm = 20*lambda_0 # Waveguide max depth
Wc =  10*lambda_0 # Width with constant depth
Wm = 35*lambda_0 # Waveguide max width. (originally infinite but we have to truncate for computational purposes).
PML_size = 4*lambda_0 # Length of the Perfectly Matched Layer (PML) that helps us truncate our computational domain.

# Location and size of the scatterer.
x_sc = 22*lambda_0 # Location in x-axis
y_sc = 8.5*lambda_0 # Location in y-axis
b = 3*lambda_0 # Size of the scatterer (radius)

## Source present in the waveguide.
frq = 73.               # Frequency in which the source emits its pulse.
omega = 2.*np.pi*frq    # Angular frequency.
k = omega/c0            # wavenumber
x_s= 4*lambda_0         # Position of source in x-axis.
y_s=  3*lambda_0        # Position of source in y-axis.
r = 5                   # Radius of source
alpha = log(10^6)/r**2  # Decay rate of the source
pulse = sqrt(alpha/pi)*exp(-alpha*((x-x_s)*(x-x_s) + (y-y_s)*(y-y_s))) # Source pulse


# We work on an infinite waveguide, which means it is bounded on top and bottom and the solution is outgoing on the left and right.
#
#           -----------------------------------------------------------------------------
#                                                               __
#                                                              /  \
#                                                              \__/
#
#           ---------------------------
#                                       \
#                                        \                   
#                                         \
#                                          \
#                                           \
#                                            --------------------------------------------



## Creating the waveguide geometry.
geo = SplineGeometry()

# Defining the points that will define the waveguide geometry.
pnts =[(0,0), #1
       (Wm,0), #2
       (Wm,Dm),  #3
       (Wm-4*lambda_0,Dm), #4
       (Wm-7*lambda_0,Dm), #5
       (0.5*(Wc+Wm-5*lambda_0), (Dc+Dm)/2), #6
       (Wc,Dc), #7
       (Wc-2*lambda_0,Dc), #8
       (0,Dc)] #9

# Appending the points to the geometry.
p1,p2,p3,p4,p5,p6,p7,p8,p9 = [geo.AppendPoint(*pnt) for pnt in pnts]

# Defining the curves that will define the waveguide geometry.
curves = [[["line",p1,p2],"top"],
          [["line",p2,p3],"right"],
          [["line",p3,p4],"bottom"],
          [["spline3",p4,p5,p6],"bottom"],
          [["spline3",p6,p7,p8],"bottom"],
          [["line",p8,p9],"bottom"],
          [["line",p9,p1],"left"]]

# Appending the curves to the geometry.
[geo.Append(c,bc=bc) for c,bc in curves]

# Adding the PMLs to the geometry.
geo.AddRectangle((-PML_size,0),(0,Dc),leftdomain=2,bc="PMLL") # Add left PML rectangle.
geo.AddRectangle((Wm,0),(Wm+PML_size,Dm),leftdomain=3,bc="PMLR") # Add right PML rectangle.
geo.AddCircle((x_sc,y_sc),b,leftdomain=0,rightdomain=1,bc="scatterer") # Add scatterer in the domain.
geo.SetMaterial(2,"PMLL")
geo.SetMaterial(3,"PMLR")

# Generating the mesh.
mesh = Mesh(geo.GenerateMesh(maxh=5))
mesh.Curve(3)

## Optional: Draw the mesh before proceeding.
# Draw(mesh)


## Letting the mesh know which are the PML layers. We then tell the mesh which labels correspond to PMLs. 
## In our case it's the PMLL and PMLR rectangles we defined earlier.
mesh.SetPML(pml.Cartesian((0,0), (Wm,Dm), 2j),"PMLL|PMLR") 

## Optional: Draw the mesh before proceeding to see the PMLs.
# #Draw(mesh); Optional meshing drawing to see our work before proceeding.

## Creating the finite element space based on the mesh we just created.
## We define the order of the finite elements, and boundary conditions (leave blank for Neummann b.c.'s)
fes = H1(mesh, order=2, complex=True, dirichlet='top|bottom|scatterer') 

u, v = fes.TnT() # Creating Test and Trial functions u, v.

Draw(pulse, mesh,'mesh') # Optional drawing to see what the source looks like.

    # Creating the weak form of the Helmholtz equation -Du - k^2 u = f
a = BilinearForm(fes, symmetric=True) # Setting a as a bilinear form
a += grad(u)*grad(v)*dx - k*k*u*v*dx
a.Assemble()

f = LinearForm(fes) # RHS is a linear form that contains the source.
f += pulse * v * dx
f.Assemble()

    # Create the grid function gfu that will contain our solution.
gfu = GridFunction(fes, name="u")

    # Solve the linear system.
gfu.vec.data = a.mat.Inverse() * f.vec

# Draw the modulus of the complex solution on the mesh.
Draw(Norm(gfu),mesh,'mesh',)

# Creating the data folder if it doesn't exist. Helps us keep the data organized and add the folder to .gitignore if needed.
if not os.path.exists(os.path.join('.', 'data')):
    os.mkdir(os.path.join('.', 'data'))

# Saving the mesh as Gmsh2 format.
meshname = "data/mesh_wg_var_depth.msh"
mesh.ngmesh.Export(meshname,"Gmsh2 Format") # Saving the mesh file. Not needed if you choose to interpolate to a regular grid later on.

# Saving the solution to a .mat file.
sol_on_mesh = ConvertSolutiononMesh(mesh,gfu) # Only keeping parts of the solution that are on mesh points and not all DOFs.

filename = "data/sol_wg_var_depth_f"+ str(frq) + ".mat"
savemat(filename,{"u":sol_on_mesh}) # Save a mat file of the Green's function on a regular grid.
