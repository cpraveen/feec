# k=1 Hodge Laplace on annular domain

The space of Harmonic forms is non-trivial and is of dimension one. We first find a basis for this in the mixed formulation. The primal formulation will converge to wrong solution.

Generate mesh

```shell
gmsh -2 annulus.geo
```

Run primal and mixed formulations

```shell
python primal.py
visit -o sol.pvd

python mixed.py
visit -o sol.pvd
```

Compare with Fig. 5.2 of Arnold
