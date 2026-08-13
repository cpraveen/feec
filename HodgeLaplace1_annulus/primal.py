'''
k = 1 Hodge Laplace problem in Gamma-shaped domain
Primal formulation, see Eq. (4.31) and (5.2). 
Boundary condition u.n = 0 is implemented weakly via penalty terms.
'''
from firedrake import *

degree = 1
gamma = 1.0e3

mesh = Mesh("annulus.msh")
print(f"Number of cells: {mesh.num_cells()}")
print(f"Number of vertices: {mesh.num_vertices()}")
print(f"Degree: {degree}")

V = VectorFunctionSpace(mesh, "Lagrange", degree)
u = TrialFunction(V)
v = TestFunction(V)

x, y = SpatialCoordinate(mesh)
n = FacetNormal(mesh)

f = as_vector([0.0,x])
a = (curl(u)*curl(v) + div(u)*div(v))*dx + gamma*dot(u,n)*dot(v,n)*ds
L = dot(f,v)*dx

u = Function(V,name="u")
solve(a == L, u)
VTKFile("sol.pvd").write(u)
