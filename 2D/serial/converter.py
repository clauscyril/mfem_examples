import meshio

# Charger le fichier .msh de Gmsh
mesh = meshio.read("data/mesh.msh", file_format="gmsh")
print(type(mesh))
# Sauvegarder en format MFEM (ASCII)
# meshio.write("data/mesh.mesh", mesh, file_format="mfem")
