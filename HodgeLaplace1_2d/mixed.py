'''
k = 1 Hodge Laplace problem in Gamma-shaped domain
Mixed formulation, see Eq. (4.32). 
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

S = FunctionSpace(mesh, "Lagrange", degree)
V = FunctionSpace(mesh, "N1curl", degree)
W = S * V
sigma,u = TrialFunctions(W)
tau,v = TestFunctions(W)

f = as_vector([1.0,0.0])
a = sigma*tau*dx - dot(u,grad(tau))*dx \
    + dot(grad(sigma),v)*dx + curl(u)*curl(v)*dx
L = dot(f,v)*dx

w = Function(W)
solve(a == L, w)
sigma,u = split(w)
u = project(u, VectorFunctionSpace(mesh, "CG", degree), name="u")
VTKFile("sol.pvd").write(u)
