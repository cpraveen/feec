'''
k = 1 Hodge Laplace problem in annular domain
Mixed formulation, see Eq. (4.32)
'''
from firedrake import *

degree = 1

mesh = Mesh("annulus.msh")
print(f"Number of cells: {mesh.num_cells()}")
print(f"Number of vertices: {mesh.num_vertices()}")
print(f"Degree: {degree}")

# Domain with hole: Harmonic space has dimension = 1. Such vectors satisfy
#    curl(q) = 0, div(q) = 0 in Omega
#    q.n = 0 on boundary
# Set q = curl(psi) and solve for psi
# q.n = 0 implies psi = const on each boundary
# Problem for psi
#    Laplace(psi) = 0, 
#    psi=0 on inner boundary
#    psi=1 on outer boundary
P = FunctionSpace(mesh, "Lagrange", degree)
p = TrialFunction(P)
q = TestFunction(P)
a = dot(grad(p),grad(q))*dx
L = Constant(0)*q*dx
bc0= DirichletBC(P, 0.0, 1) # inner boundary
bc1= DirichletBC(P, 1.0, 2) # outer boundary
psi = Function(P)
solve(a == L, psi, bcs=[bc0,bc1])

# Basis for Harmonic space
q = curl(psi)

S = FunctionSpace(mesh, "Lagrange", degree)
V = FunctionSpace(mesh, "N1curl", degree)
R = FunctionSpace(mesh, "R", 0)
W = S * V * R
sigma,u,r = TrialFunctions(W)
tau,v,s   = TestFunctions(W)

# All bc are natural

x, y = SpatialCoordinate(mesh)
f = as_vector([0.0,x])
a = sigma*tau*dx - dot(u,grad(tau))*dx \
    + dot(grad(sigma),v)*dx + curl(u)*curl(v)*dx \
    + r*dot(q,v)*dx + s*dot(u,q)*dx
L = dot(f,v)*dx

w = Function(W)
solve(a == L, w)
sigma,u,r = split(w)
u = project(u, VectorFunctionSpace(mesh, "CG", degree), name="u")
VTKFile("sol.pvd").write(u)
