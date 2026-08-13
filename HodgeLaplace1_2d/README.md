# k=1 Hodge Laplace problem in 2d

Make the mesh

```shell
python mesh.py
```

See the mesh

```shell
gmsh gamma.msh
```

Run primal formulation

```shell
python primal.py
visit -o sol.pvd
```

and plot vectors of `u`.

Run mixed formulation

```shell
python mixed.py
visit -o sol.pvd
```

and plot vectors of `u`.

<img src="output/mixed.png" alt="Mixed"  width="45%"/>
<img src="output/primal.png" alt="Primal"  width="45%"/>

Compare with Fig. 5.1 of Arnold.
