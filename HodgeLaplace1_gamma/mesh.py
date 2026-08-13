import gmsh

def generate_gamma_mesh(mesh_size=0.05, filename="gamma.msh"):
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)  # Suppress console output
    gmsh.model.add("gamma")

    # Define vertices (Points)
    lc = mesh_size
    p1 = gmsh.model.geo.addPoint(0.0, 0.5, 0.0, lc)
    p2 = gmsh.model.geo.addPoint(0.5, 0.5, 0.0, lc)
    p3 = gmsh.model.geo.addPoint(0.5, 0.0, 0.0, lc)
    p4 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0, lc)
    p5 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0, lc)
    p6 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0, lc)

    # Define boundary segments (Lines)
    l1 = gmsh.model.geo.addLine(p1, p2) # 2
    l2 = gmsh.model.geo.addLine(p2, p3) # 1
    l3 = gmsh.model.geo.addLine(p3, p4) # 2
    l4 = gmsh.model.geo.addLine(p4, p5) # 1
    l5 = gmsh.model.geo.addLine(p5, p6) # 2
    l6 = gmsh.model.geo.addLine(p6, p1) # 1

    # Create curve loop and surface
    loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4, l5, l6])
    surface = gmsh.model.geo.addPlaneSurface([loop])

    gmsh.model.geo.synchronize()

    # Assign physical IDs for boundaries
    gmsh.model.addPhysicalGroup(1, [l2, l4, l6], tag=1, name="xnormal")
    gmsh.model.addPhysicalGroup(1, [l1, l3, l5], tag=2, name="ynormal")

    # Assign physical ID for 2D domain surface
    gmsh.model.addPhysicalGroup(2, [surface], tag=100, name="domain")

    # Generate 2D triangular mesh
    gmsh.model.mesh.generate(2)
    gmsh.write(filename)
    gmsh.finalize()


if __name__ == "__main__":
    # Generate mesh with characteristic element length h = 0.05
    generate_gamma_mesh(mesh_size=0.01)
    print(f"Mesh created successfully.")
