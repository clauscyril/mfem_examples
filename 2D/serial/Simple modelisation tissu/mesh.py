import meshpy.triangle as triangle
import numpy as np

def circle_points(n, radius=1.0):
    """Génère des points sur le contour d'un disque."""
    theta = np.linspace(0, 2*np.pi, n, endpoint=False)
    return [(radius*np.cos(t), radius*np.sin(t)) for t in theta]

# 1. Génération du maillage avec MeshPy
boundary_points = circle_points(40)
facets = [(i, (i+1) % len(boundary_points)) for i in range(len(boundary_points))]

info = triangle.MeshInfo()
info.set_points(boundary_points)
info.set_facets(facets)

mesh = triangle.build(info, max_volume=0.005)

# 2. Exportation au format MFEM .mesh
def export_to_mfem(filename, mesh):
    with open(filename, "w") as f:
        f.write("MFEM mesh v1.0\n\n")
        f.write("dimension\n2\n\n")

        # Éléments (triangles)
        f.write(f"elements\n{len(mesh.elements)}\n")
        for i, e in enumerate(mesh.elements):
            f.write(f"1 2 {e[0]} {e[1]} {e[2]}\n")

        # Frontière (segments du bord)
        f.write(f"\nboundary\n{len(mesh.facets)}\n")
        for i, (nodes, _) in enumerate(mesh.facets):
            f.write(f"1 1 {nodes[0][0]} {nodes[0][1]}\n")

        # Sommets
        f.write(f"\nvertices\n{len(mesh.points)}\n2\n")
        for p in mesh.points:
            f.write(f"{p[0]} {p[1]}\n")

export_to_mfem("disk.mesh", mesh)
print("Maillage exporté sous 'disk.mesh'.")
