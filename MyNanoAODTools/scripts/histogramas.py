import uproot
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys


#Directorio con archivos root como argumento
if len(sys.argv) < 2:
    sys.exit(1)
    
input_directory = sys.argv[1]

# Ruta a los archivos .root
root_files = glob.glob(os.path.join(input_directory,"*Skim.root"))

# Nombre del TTree
tree_name = "Events"

# Lista de branches a procesar
branch_list = ["A_Zmass", "A_Dr_Z", "A_Sum_pt", "A_Sum_mass",
               "B_Zmass", "B_Dr_Z", "B_Sum_pt", "B_Sum_mass",
               "C_Zmass", "C_Dr_Z", "C_Sum_pt", "C_Sum_mass",
               "D_Zmass", "D_Dr_Z", "D_Sum_pt", "D_Sum_mass"]  # <-- Agrega aqui los branches

# Obtener directorio base de los archivos
output_dir = os.path.dirname(root_files[0]) if root_files else "."

# Procesar cada branch
for branch_name in branch_list:
    all_data = []

    for filename in root_files:
        with uproot.open(filename) as file:
            tree = file[tree_name]

            # Verificar si el branch existe en ese archivo
            if branch_name in tree.keys():
                branch_data = tree[branch_name].array(library="np")
                all_data.append(branch_data)
            else:
                print(f"[WARN] '{branch_name}' no encontrado en {filename}")

    if not all_data:
        print(f"[SKIP] No se encontraron datos para '{branch_name}' en ningun archivo.")
        continue

    # Combinar todos los datos en uno solo
    combined_data = np.concatenate(all_data)

    # Crear histograma
    plt.figure(figsize=(10,6))
    plt.hist(combined_data, bins=100, histtype='stepfilled', alpha=0.7)
    plt.xlabel(branch_name)
    plt.ylabel("Entries")
    plt.title(f"Histograma combinado - Branch '{branch_name}'")
    plt.grid(True)

    # Guardar el histograma
    output_path = os.path.join(output_dir, f"hist_{branch_name}.png")
    plt.savefig(output_path)
    plt.close()  # Importante para evitar sobrecarga de memoria al hacer multiples figuras

    print(f"[OK] Histograma guardado en: {output_path}")
