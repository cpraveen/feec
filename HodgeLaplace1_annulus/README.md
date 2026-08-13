# k=1 Hodge Laplace on annular domain

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
