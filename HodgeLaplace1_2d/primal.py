'''
k = 1 Hodge Laplace problem in Gamma-shaped domain
Primal formulation, see Eq. (4.31) and (5.2). 
'''
from firedrake import *

filename = "gamma.msh"
degree = 1

mesh = Mesh(filename)

print(f"Number of cells: {mesh.num_cells()}")
print(f"Number of vertices: {mesh.num_vertices()}")
print(f"Degree: {degree}")

V = VectorFunctionSpace(mesh, "Lagrange", degree)
u = TrialFunction(V)
v = TestFunction(V)

# u.n = 0
# tag=1 are vertical faces
# tag=2 are horizontal faces
bcx = DirichletBC(V.sub(0), 0.0, 1)
bcy = DirichletBC(V.sub(1), 0.0, 2)
bcs = [bcx, bcy]

f = as_vector([1.0,0.0])
a = (curl(u)*curl(v) + div(u)*div(v))*dx
L = dot(f,v)*dx

u = Function(V,name="u")
solve(a == L, u, bcs=bcs)
VTKFile("sol.pvd").write(u)
