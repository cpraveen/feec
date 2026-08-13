'''
k = 1 Hodge Laplace problem in Gamma-shaped domain
Primal formulation, see Eq. (4.31) and (5.2). 
p=0 so we do not need to solve for this.
'''
from firedrake import *
import matplotlib.pyplot as plt

filename = "gamma.msh"
degree = 1

mesh = Mesh(filename)

print(f"Number of cells: {mesh.num_cells()}")
print(f"Number of vertices: {mesh.num_vertices()}")
print(f"Degree: {degree}")

V = VectorFunctionSpace(mesh, "Lagrange", degree)
u = TrialFunction(V)
v = TestFunction(V)

bcx = DirichletBC(V.sub(0), 0.0, 1)
bcy = DirichletBC(V.sub(1), 0.0, 2)
bcs = [bcx, bcy]

f = as_vector([1.0,0.0])
a = (curl(u)*curl(v) + div(u)*div(v))*dx
L = dot(f,v)*dx

u = Function(V,name="u")
solve(a == L, u, bcs=bcs)
VTKFile("sol.pvd").write(u)
