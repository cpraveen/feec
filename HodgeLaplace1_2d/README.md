# k=1 Hodge Laplace problem in 2d

Solve in simply connected domain

```math
\textrm{curl}(\textrm{curl}(u)) - \textrm{grad}(\textrm{div}(u)) = f
```

with boundary condition

```math
u \cdot n = 0, \qquad (\textrm{curl } u) \times n = 0
```

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

Compare with Fig. 5.1 of Arnold. Left figures shows mixed formulation and right shows the primal forumation. The primal formulation converges to wrong solution.
