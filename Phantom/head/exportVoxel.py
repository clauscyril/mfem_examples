import numpy as np
import scipy.io

# Charger le fichier .mat
mat_data = scipy.io.loadmat("tete.mat")

# Extraire le volume (remplace 'voxel' par le vrai nom dans ton .mat)
volume = mat_data["data"]  # Doit être un ndarray 3D de dtype uint8

print(mat_data["data"])

# Vérifier la forme
ZDIM, YDIM, XDIM = volume.shape

# Créer l’en-tête INR
header = f"""#INRIMAGE-4#{{
XDIM={XDIM}
YDIM={YDIM}
ZDIM={ZDIM}
VDIM=1
TYPE=unsigned char
PIXSIZE=8 bits
SCALE=2**0
CPU=decm
VX=1.0
VY=1.0
VZ=1.0
##}}
"""

# Sauvegarder l’image INR
with open("imTete.inr", "wb") as f:
    f.write(header.encode("binary"))
    f.write(volume.tobytes(order="C"))  # Ordre mémoire : X → Y → Z
