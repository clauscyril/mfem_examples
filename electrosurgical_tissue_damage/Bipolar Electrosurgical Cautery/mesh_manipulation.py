import os

# Chemin du fichier source
path = os.path.dirname(os.path.abspath(__file__))
source_file = os.path.join(path, "meshOnshape_test.msh")

# Lecture du fichier .msh
with open(source_file, "r", encoding="utf-8") as f:
    texte = f.read()

# Localisation des sections
idx_nodes_start = texte.find("$Nodes")
idx_elements_start = texte.find("$Elements")

header = texte[:idx_nodes_start]
nodes_section = texte[idx_nodes_start:idx_elements_start]
elements_section = texte[idx_elements_start:]

# Extraction des lignes utiles
Nodes = nodes_section.splitlines()[2:-1]      # Ignorer $Nodes et $EndNodes
Elements = elements_section.splitlines()[2:-1]  # Ignorer $Elements et $EndElements

# Régions physiques à supprimer
regions_to_delete = {"85", "86", "90", "91"}

# Filtrage des éléments et collecte des nœuds à supprimer
nodes_to_delete = set()
filtered_elements = []

for line in Elements:
    parts = line.strip().split()
    if parts[3] in regions_to_delete or parts[4] in regions_to_delete:
        nodes_to_delete.update(parts[5:])
    else:
        filtered_elements.append(parts)

# Renumérotation des nœuds conservés
node_map = {}
new_nodes = []
counter = 1

for line in Nodes:
    parts = line.strip().split()
    node_id = parts[0]
    if node_id not in nodes_to_delete:
        node_map[node_id] = str(counter)
        new_nodes.append([str(counter)] + parts[1:])
        counter += 1

# Mise à jour des éléments avec les nouveaux IDs de nœuds
new_elements = []
for i, parts in enumerate(filtered_elements, start=1):
    updated_nodes = [node_map[n] for n in parts[5:]]
    new_elements.append([str(i)] + parts[1:5] + updated_nodes)

# Reconstruction du fichier .msh
new_nodes_text = "\n".join(" ".join(n) for n in new_nodes)
new_elements_text = "\n".join(" ".join(e) for e in new_elements)

new_file_content = (
    f"{header}$Nodes\n{new_nodes_text}\n$EndNodes\n"
    f"$Elements\n{new_elements_text}\n$EndElements"
)

# Écriture du nouveau fichier
output_file = os.path.join(path, "new_mesh2.msh")
with open(output_file, "w", encoding="utf-8") as f:
    f.write(new_file_content)

print(f"Nouveau fichier sauvegardé sous : {output_file}")
